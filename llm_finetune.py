### 对文献进行训练
### 
import wandb
from transformers import (
    LlamaForCausalLM, LlamaTokenizer,
    AutoModel, AutoTokenizer, AutoModelForCausalLM,
    BloomForCausalLM, BloomTokenizerFast
)
from save import SavePeftModelCallback

from peft import (
    prepare_model_for_int8_training,
    AdaLoraConfig,
    PrefixTuningConfig,
    PromptEncoderConfig,
    PromptTuningInit,
    LoraConfig,
    get_peft_model,
    get_peft_model_state_dict,
)
import copy
import os
import argparse
from collections import namedtuple

ModelClass = namedtuple("ModelClass", ('tokenizer', 'model'))

_MODEL_CLASSES = {
    "llama": ModelClass(**{
        "tokenizer": LlamaTokenizer,
        "model": LlamaForCausalLM,
    }),
    "bloom": ModelClass(**{
        "tokenizer": BloomTokenizerFast,
        "model": BloomForCausalLM,
    }),
    "Auto": ModelClass(**{
        "tokenizer": AutoTokenizer,
        "model": AutoModel,
    }),

}

_PEFT_CLASSES = {
    "lora": LoraConfig,
    "adalora": AdaLoraConfig,
}
device_map = "auto"
# device_map，是将模型切分为很多块，使得
world_size = int(os.environ().get("WORLD_SIZE", 1))
ddp = world_size != 1

import transformers

# add the custom dataset
DATA_PATH = {
    "neurology": "pubmed_data_neurology.json",
}

IGNORE_INDEX = -100

LOAD_8BIT = True

PROMPT_DICT = {
    "prompt_input": (
        "Below is an instruction that describes a task, paired with an input that provides further context. "
        "Write a response that appropriately completes the request.\n\n"
        "### Instruction:\n{instruction}\n\n### Input:\n{input}\n\n### Response:"
    ),
    "prompt_no_input": (
        "Below is an instruction that describes a task. "
        "Write a response that appropriately completes the request.\n\n"
        "### Instruction:\n{instruction}\n\n### Response:"
    ),
    "prompt_multirun_input": (
        "Below is an multi-round dialogue between human and assistant. "
        "Write a response as an assistant that appropriately completes the human request in each round by incorporating previous context.\n\n"
        "{instruction}{output}"
    ),
}

def generate_prompt(data_point):
    # a nasty solution just for now
    prompt_ = PROMPT_DICT['prompt_input'] if data_point["input"] else PROMPT_DICT['prompt_no_input']
    return prompt_.format_map(data_point)

from datasets import load_dataset, concatenate_datasets, DatasetDict

def get_data_model(args):
    # data, model, tokenizer，通过一个函数就获取了
    def get_model_class(model_type):
        if model_type not in ['bloom', 'llama', 'chatglm', 'moss']:
            model_type = "Auto"
            
        
        return _MODEL_CLASSES[model_type] # tokenizer, model
    
    def get_peft_class(peft_type):

        return _PEFT_CLASSES[peft_type] # tokenizer, model

    data = DatasetDict()

    if len(args.data) == 1 and not  args.data[0].endswith(".json"):
        data_file_path = DATA_PATH.get(args.data[0], None)
        assert data_file_path, "Error: Wrong type of data."
        data = load_dataset("json", data_files=data_file_path)
    else:
        merge_data = concatenate_datasets([load_dataset("json", data_files=fname)["train"] for fname in args.data])
        data = DatasetDict({"train":merge_data})

    print(data)

    model_class = get_model_class(args.model_type)
    peft_class = get_peft_class(args.peft_type)

    if args.model_type in ["chatglm"]:
        model = model_class.model.from_pretrained(args.model_name_or_path,
                                                trust_remote_code=True,
                                                device_map=device_map)
        tokenizer = model_class.tokenizer.from_pretrained(args.model_name_or_path, trust_remote_code=True)
    elif args.model_type in ["moss"]:
        model = model_class.model.from_pretrained(args.model_name_or_path,
                                                trust_remote_code=True,
                                                load_in_8bit=LOAD_8BIT,
                                                device_map=device_map)
        tokenizer = model_class.tokenizer.from_pretrained(args.model_name_or_path, trust_remote_code=True)
    else:
        model = model_class.model.from_pretrained(args.model_name_or_path,
                                                load_in_8bit=LOAD_8BIT,
                                                device_map=device_map)
        tokenizer = model_class.tokenizer.from_pretrained(args.model_name_or_path) # default add_eos_token=False

    # llama has no pad_id, maybe copy the stanford_alpaca's handling ?
    if args.model_type in ['llama', 'moss']:
        tokenizer.pad_token_id = 0 # unk_id in llama. we want this to be different from the eos token

    if LOAD_8BIT:
        model = prepare_model_for_int8_training(model)

    if args.peft_type=='lora':
        # 据说lora是要为每个模型专门写一个
        config = peft_class(
            init_r=args.lora_r,
            lora_alpha=args.lora_alpha,
            target_modules=args.lora_target_modules,
            lora_dropout=args.lora_dropout,
            bias="none",
            task_type="CAUSAL_LM",
        )
    elif args.peft_type=="adalora":
        config = peft_class(
            init_r=args.adalora_init_r,
            r=args.lora_r,
            beta1=0.85,
            beta2=0.85,
            tinit=args.adalora_tinit,
            tfinal=args.adalora_tfinal,
            deltaT=args.adalora_delta_t,
            lora_alpha=args.lora_alpha,
            lora_dropout=args.lora_dropout,
            target_modules=args.lora_target_modules,
            task_type="CAUSAL_LM",
            inference_mode=False,
        )
    else:
        assert args.peft_type, "Error: Wrong type of peft."

    model = get_peft_model(model, config)

    # the size of trainable parameters for lora modules
    model.print_trainable_parameters()

    return data, model, tokenizer


def train(args):
    # 1. load data & model_class
    data, model, tokenizer = get_data_model(args)
    
    if "chatglm" in args.model_type:
        def prompt_tokenize(prompt):
            input_ids = tokenizer.encode(prompt,
                                         truncation=True,
                                         max_length=args.cutoff_len,
                                         padding=False,
                                         )# 为什么不设置padding=True
            return {
                "input_ids": input_ids,
                "labels": copy.deepcopy(input_ids),
            }
            # 这里的labels不用管device
            
        def completion_tokenize(completion):
            input_ids = tokenizer.encode(completion, max_length=args.cutoff_len)
            # 这里需要把不满足长度的，填充到最长
            return {
                "input_ids": input_ids,
                "labels": copy.deepcopy(input_ids)
            }

    elif "moss" in args.model_type:
        def tokenize(prompt):
            result = tokenizer(
                prompt,
                truncation=True,
                max_length=args.cutoff_len,
                # padding="max_length",
            )
            return {
                "input_ids": result["input_ids"],
                "labels": copy.deepcopy(result["input_ids"]),
                "attention_mask": result["attention_mask"],
            }
    else:
        def tokenize(prompt):
            result = tokenizer(prompt,
                               truncation=True,
                               max_length=args.cutoff_len,
                            #    padding="max_length",
                               padding=False,
                            )
            return {
                "input_ids": result["input_ids"],
                "attention_mask": result["attention_mask"],
                "labels": copy.deepcopy(result["input_ids"])
            }

    def generate_and_tokenize_prompt(data_point):
        prompt_no_resp = generate_prompt(data_point)
        if "chatglm" in args.model_type:
            tokenized_result = prompt_tokenize(prompt_no_resp)
        else:
            tokenized_result = tokenize(prompt_no_resp)

        source_len = len(tokenized_result['input_ids'])
        prompt_with_response = prompt_no_resp + " " + data_point["output"]
        # if "llama" in args.model_type:
        prompt_with_response += " " + tokenizer.eos_token
        if "chatglm" in args.model_type:
            tokenized_with_response = completion_tokenize(prompt_with_response)
        else:
            tokenized_with_response = tokenize(prompt_with_response)
        tokenized_with_response["labels"] = [IGNORE_INDEX] * source_len + tokenized_with_response["labels"][source_len:]

        return tokenized_with_response

    model_name = args.model_name_or_path.split( '/')[-1]
    data_name = "+".join([d.split("/")[-1].strip(".json") for d in args.data])
    lr_str = str(args.learning_rate)
    output_dir = f"saved_models/{model_name}_{data_name}_{lr_str}/{args.peft_type}"

    wandb.init(
        project = "Alpaca-CoT",
        config = args,
        name = f"{model_name}_{data_name}_{lr_str}"
    )

    # 2. split dataset
    if args.val_set_size > 0:
        train_val = data["train"].train_test_split(
            test_size=args.val_set_size, shuffle=True, seed=42
        )
        train_data = train_val["train"].shuffle().map(generate_and_tokenize_prompt)
        val_data = train_val["test"].shuffle().map(generate_and_tokenize_prompt)
    else:
        train_data = data["train"].shuffle().map(generate_and_tokenize_prompt)
        val_data = None

    # 3. train
    total_batch_size = args.per_gpu_train_batch_size * args.gradient_accumulation_steps * (world_size if ddp else 1)
    total_optim_steps = train_data.num_rows // total_batch_size
    saving_step = int(total_optim_steps/10)
    warmup_steps = int(total_optim_steps/10)

    print("***** Running training *****")
    print(f"  Num Epochs = {args.epochs}", )
    print(f"  Instantaneous batch size per GPU = {args.per_gpu_train_batch_size}")
    print(f"  Gradient Accumulation steps = {args.gradient_accumulation_steps}")
    print(f"  Total train batch size (w. parallel, distributed & accumulation) = {total_batch_size}")
    print(f"  Total optimization steps = {total_optim_steps}")
    print(f"  Saving steps = {saving_step}")

    trainer = transformers.Trainer(
        model=model,
        train_dataset=train_data,
        eval_dataset=val_data,
        args=transformers.TrainingArguments(
            per_device_train_batch_size=args.per_gpu_train_batch_size,
            gradient_accumulation_steps=args.gradient_accumulation_steps,
            warmup_steps=warmup_steps,
            num_train_epochs=args.epochs,
            learning_rate=args.learning_rate,
            fp16=True,
            logging_steps=20,
            evaluation_strategy="steps" if args.val_set_size > 0 else "no",
            save_strategy="steps",
            eval_steps=saving_step if args.val_set_size > 0 else None,
            save_steps=saving_step,
            output_dir=output_dir,
            save_total_limit=11,
            load_best_model_at_end=True if args.val_set_size > 0 else False,
            ddp_find_unused_parameters=False if ddp else None,
        ),
        data_collator=transformers.DataCollatorForSeq2Seq(tokenizer, return_tensors="pt", padding=True),
        callbacks=[SavePeftModelCallback],
    )
    model.config.use_cache = False
    if torch.__version__ >= "2" and sys.platform != "win32":
        model = torch.compile(model)

    trainer.train(resume_from_checkpoint=args.resume_from_checkpoint)

    model.save_pretrained(output_dir)

    print("\n If there's a warning about missing keys above, please disregard :)")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--size', type=str, help='the size of llama model')
    parser.add_argument('--data', type=str, nargs="*", help='the data used for instructing tuning')
    parser.add_argument('--local_rank', default=-1, type=int, help='node rank for distributed training')
    parser.add_argument('--model_type', default="llama", choices=['llama', 'chatglm', 'bloom', 'moss'])
    parser.add_argument('--model_name_or_path', default="decapoda-research/llama-7b-hf", type=str)
    parser.add_argument('--per_gpu_train_batch_size', default=4, type=int, help='Batch size per GPU/CPU for training.')
    parser.add_argument('--gradient_accumulation_steps', default=32, type=int)
    parser.add_argument('--epochs', default=3, type=int)
    parser.add_argument('--learning_rate', default=3e-4, type=float)
    parser.add_argument('--cutoff_len', default=512, type=int)
    #PEFT arguments
    parser.add_argument('--peft_type', default="lora", choices=['lora', 'adalora', 'prompt','p_tuning','prefix'])
    parser.add_argument('--lora_r', default=8, type=int)
    parser.add_argument('--lora_alpha', default=16, type=int)
    parser.add_argument('--lora_dropout', default=0.05, type=float)
    parser.add_argument('--val_set_size', default=2000, type=int)
    parser.add_argument('--lora_target_modules', nargs='+',
                        help="the module to be injected, e.g. q_proj/v_proj/k_proj/o_proj for llama, query_key_value for bloom&GLM",
                        default=["q_proj", "v_proj"])
    parser.add_argument('--adalora_init_r', default=12, type=int)
    parser.add_argument("--adalora_tinit", type=int, default=200, help="number of warmup steps for AdaLoRA wherein no pruning is performed")
    parser.add_argument("--adalora_tfinal", type=int, default=1000, help=" fix the resulting budget distribution and fine-tune the model for tfinal steps when using AdaLoRA ")
    parser.add_argument("--adalora_delta_t", type=int, default=10, help="interval of steps for AdaLoRA to update rank")
    parser.add_argument('--num_virtual_tokens', default=20, type=int)
    parser.add_argument('--prompt_encoder_hidden_size', default=128, type=int)
    parser.add_argument('--resume_from_checkpoint', nargs='?', default=None, const=True, help='resume from the specified or the latest checkpoint, e.g. `--resume_from_checkpoint [path]` or `--resume_from_checkpoint`')

    args, _ = parser.parse_known_args()
    print(args)

    train(args)