# 翻译

import json

#from Baidu_Text_transAPI import translate
from ChatGPT_Translate import translate

end_year=2021
start_year=2000
term="psychiatry"

# baidu api 5XBoy86ucIaqjsHQfSZc
from tqdm import tqdm
for year in range(end_year, start_year, -1):
    source_file = f"{term}LLM/papers_{term}_{year}.json"
    with open(source_file, "r") as f:
        papers = json.load(f)
    translate_json_list = []
    for index, paper in enumerate(tqdm(papers)):
        if index % 10 == 0:
            translate_json = {}
            translate_json["Title"] = translate(paper["Title"])
            translate_json["Abstract"] = translate(paper["Abstract"])
            translate_json_list.append(translate_json)
            print(translate_json)
        
    # Save the results to a json file
    with open("Zh"+source_file, "w") as zhf:
        json.dump(translate_json_list, zhf)
