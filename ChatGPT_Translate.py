
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