import time
import json
from Bio import Entrez

# 设置你的email和可选的API key
Entrez.email = "cui_yaxin@outlook.com"
Entrez.api_key = "b18fac6e0128aa51c1b86a6b047071f37a08"

# 搜索关键词
term = "psychiatry"

# 每次请求的记录数
batch_size = 200

# 总的记录数
total_records = 100000

# 保存数据的文件
output_file = "pubmed_data.json"

# 初始化起始位置
start = 0

# 保存所有的数据
all_data = []

while start < total_records:
    try:
        # 使用esearch函数搜索论文
        search_results = Entrez.read(Entrez.esearch(db="pubmed", term=term, retmax=batch_size, retstart=start, sorted="pub date"))

        # 获取搜索结果的PMID列表
        pmid_list = search_results["IdList"]
        
        # 使用efetch函数获取PMID对应的详细信息
        fetch_results = Entrez.efetch(db="pubmed", id=pmid_list, rettype="xml")

        # 解析返回的XML数据
        papers = Entrez.read(fetch_results)

        # 处理每篇论文的数据
        for paper in papers["PubmedArticle"]:
            # 获取题目
            print("paper")
            print(paper)
            print("--"*100)
            title = paper["MedlineCitation"]["Article"]["ArticleTitle"]
            
            # 检查是否有摘要
            if "Abstract" in paper["MedlineCitation"]["Article"]:
                abstract = paper["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
                # 假设Abstract的第一部分就是Introduction
                introduction = abstract[0] if abstract else ""
            else:
                abstract = ""
                introduction = ""
            # 保存数据
            all_data.append({
                "title": title,
                "abstract": abstract,
                "introduction": introduction,
            })

        # 更新起始位置
        start += batch_size
        if start % 500 ==0:
            print(all_data[-1])
        # 休眠一下，避免请求过于频繁
        time.sleep(1)

    except Exception as e:
        print("Error:", e)
        print("Retrying...")

        # 如果出错，休眠一下再重试
        time.sleep(10)

# 保存所有的数据到文件
with open(output_file, "w") as f:
    json.dump(all_data, f)

