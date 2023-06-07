
import openai


openai.api_key = "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"
import time

def translate(text="This is an English Sentence"):
    prompt = "请将以下句子翻译为中文：\n"
    prompt = prompt + text
    # Generate a response
    completion = openai.ChatCompletion.create(
        model="gpt-3.5-turbo",
        messages=[
            {"role": "user", "content": prompt}
        ]
    )
    response = completion.choices[0].message.content
    time.sleep(0.2)
    return(response.strip())