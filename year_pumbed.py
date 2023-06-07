import json
import time
from Bio import Entrez
from tqdm import tqdm

# 设置你的email和可选的API key
Entrez.email = "cui_yaxin@outlook.com"
Entrez.api_key = "b18fac6e0128aa51c1b86a6b047071f37a08"

# 每次请求的记录数
batch_size = 400

from math import sqrt, log2

# Start and end years for the search

start_year = 2000
end_year = 2023
term="kidney"
# 初始化起始位置
import os
os.mkdir(term)

id_list_set = set()
for year in range(end_year, start_year, -1):
    # Perform the search
    year_all_data = []
    # 总的记录数
    start = 0

    num = 0
    total_records = 10000*int(sqrt(start_year-1999))
    while start < total_records:
        print(start)
        try:
            query = f"{term}[Title/Abstract] AND {year}[pdat]"
            handle = Entrez.esearch(db="pubmed", term=query, retmax=batch_size, retstart=start)
            record = Entrez.read(handle)
            handle.close()

            id_list = record["IdList"]

            if not id_list:
                print(f"{year} 存入完成")
                print("--"*100)
                break
            # Fetch the details for each id
            results = []
            handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
            print("load id list success")
            # To prevent from overloading the server with requests, a delay is used.
            time.sleep(0.15)

            papers = Entrez.read(handle)
            print("read id list success")

            handle.close()

            # Extract the necessary information
            for paper in tqdm(papers["PubmedArticle"]):
                # 获取题目
                num+=1
                title = paper["MedlineCitation"]["Article"]["ArticleTitle"]
                # 检查是否有摘要
                if "MedlineCitation" not in paper.keys():
                    print("no MedlineCitation")
                    continue
                
                if "Abstract" in paper["MedlineCitation"]["Article"]:
                    AbstractText = paper["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
                    abstract = '\n'.join(AbstractText)
                else:
                    # 如果连摘要都没有，没必要训练
                    continue
                
                year_all_data.append({
                    "Title": title,
                    "Abstract": abstract,
                })
                
            # 更新起始位置
            start += batch_size
            print(num, '/', start, '/', len(year_all_data))
        except Exception as e:
            print("Error:", e)
            print("Retrying...")

            # 如果出错，休眠一下再重试
            time.sleep(10)

            
    # Save the results to a json file
    with open(f"{term}/papers_{term}_{year}.json", "w") as f:
        json.dump(year_all_data, f)
