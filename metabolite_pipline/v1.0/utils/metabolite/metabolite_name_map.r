#!/home/genesky/software/r/3.5.1/bin/Rscript
library(docopt)
"Usage: metabolite_name_map.r -i <file>  --compound_db <file> --syn_nms <file>  -o <file> 
Options:
   -i, --input <file>         代谢物名称列表文件，只看第一列，没有表头
   -c, --compound_db <file>   从metaboanalyst下载的数据库 metaboanalyst.compound_db.rds
   -s, --syn_nms <file>       从metaboanalyst下载的数据库 metaboanalyst.syn_nms.rds
   -o, --output <file>        输出文件， 例如： ./metabolite_name_map.csv " -> doc

opts                     <- docopt(doc, version = 'Program : 代谢物样本水平归一化 v1.0 \n          甘斌 129\n')
input        <- opts$input
compound_db  <- opts$compound_db
syn_nms      <- opts$syn_nms
output       <- opts$output
message("start metabolite_name_map.r")
# 以下代码框架、数据库全部来自于网站： https://www.metaboanalyst.ca/MetaboAnalyst

# input = "./metabolite_database/test.metabolite.txt"
# compound_db = "./metabolite_database/metaboanalyst.compound_db.rds"
# syn_nms = "./metabolite_database/metaboanalyst.syn_nms.rds"
message("loading")
# 读入文件、加载rds数据
data_raw <- read.table(input, head = F, check.names = F, stringsAsFactors = F, sep = "\t", quote = "")

# 注：以下两个数据库行数完全相同，且一一对应
cmpd.db <- readRDS(compound_db)
syn.db <- readRDS(syn_nms)  # 包含了代谢物的各种同义名称

# 开始注释
metabolite_names = data_raw[, 1]

# （1）初始化数据，用于记录输入代谢物与数据库的匹配情况
message("initial values")
hit.inx <- vector(mode='numeric', length=length(metabolite_names)); # record hit index, initial 0
names(hit.inx) <- metabolite_names;
match.values <- vector(mode='character', length=length(metabolite_names)); # the best matched values (hit names), initial ""
match.state <- vector(mode='numeric', length=length(metabolite_names));  # match status - 0, no match; 1, exact match; initial 0 
  
# （2）与数据库中的name完整匹配
message("exact match")
hit.inx <- match(tolower(metabolite_names), tolower(cmpd.db$name));
match.values <- cmpd.db$name[hit.inx];
match.state[!is.na(hit.inx)] <- 1;

# （3）对上一步匹配失败的代谢物，从同义数据一个一个地匹配
message("synanyms match")
syns.list <-  syn.db$syns.list;  
todo.inx <-which(is.na(hit.inx));  # 需要继续搜索的代谢物行号
if(length(todo.inx) > 0)
{
   for(i in 1:length(syns.list))
   {
        syns <-  syns.list[[i]];
        hitInx <- match(tolower(metabolite_names[todo.inx]), tolower(syns));
        
        # 如果能够匹配上，把行号记录下来
        hitPos <- which(!is.na(hitInx));
        if(length(hitPos)>0)
        {
          # record matched ones
          orig.inx<-todo.inx[hitPos];
          hit.inx[orig.inx] <- i;                  
          # match.values[orig.inx] <- syns[hitInx[hitPos]];  # show matched synnames
          match.values[orig.inx] <- cmpd.db$name[i];    # show common name
          match.state[orig.inx] <- 1;
          
          # update unmatched list
          todo.inx <- todo.inx[is.na(hitInx)];
        }
        if(length(todo.inx) == 0) break;
   }
}
 
# （4）结果整理输出
message("output")
result <- matrix("", nrow=length(metabolite_names), ncol=9);
colnames(result) <- c("Query", "Match", "HMDB", "PubChem", "ChEBI", "KEGG", "METLIN", "SMILES", "Comment");
  
for (i in 1:length(metabolite_names))
{
   #  if(match.state[i] == 1){
   #    pre.style<-"";
   #    post.style="";
   #  }else{ # no matches
   #    pre.style<-no.prestyle;
   #    post.style<-no.poststyle;
   #  }
    hit <-cmpd.db[hit.inx[i], ,drop=F];  # 去除当前代谢物的注释信息
   #  html.res[i, ]<-c(paste(pre.style, metabolite_names[i], post.style, sep=""),
   #                   paste(ifelse(match.state[i]==0, "", match.values[i]), sep=""),
   #                   paste(ifelse(match.state[i]==0 || is.na(hit$hmdb_id) || hit$hmdb_id=="" || hit$hmdb_id=="NA","-", paste("<a href=http://www.hmdb.ca/metabolites/", hit$hmdb_id, " target='_blank'>",hit$hmdb_id,"</a>", sep="")),  sep=""),
   #                   paste(ifelse(match.state[i]==0 || is.na(hit$pubchem_id) || hit$pubchem_id=="" || hit$pubchem_id=="NA", "-", paste("<a href=http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=", hit$pubchem_id," target='_blank'>", hit$pubchem_id,"</a>", sep="")), sep=""),
   #                   paste(ifelse(match.state[i]==0 || is.na(hit$chebi_id) || hit$chebi_id==""|| hit$chebi_id=="NA","-", paste("<a href=http://www.ebi.ac.uk/chebi/searchId.do?chebiId=", hit$chebi_id, " target='_blank'>",hit$chebi_id,"</a>", sep="")), sep=""),
   #                   paste(ifelse(match.state[i]==0 || is.na(hit$kegg_id) || hit$kegg_id==""|| hit$kegg_id=="NA","-",paste("<a href=http://www.genome.jp/dbget-bin/www_bget?", hit$kegg_id, " target='_blank'>", hit$kegg_id,"</a>", sep="")), sep=""),
   #                   paste(ifelse(match.state[i]==0 || is.na(hit$metlin_id) || hit$metlin_id==""|| hit$metlin_id=="NA","-",paste("<a href=http://metlin.scripps.edu/metabo_info.php?molid=", hit$metlin_id," target='_blank'>",hit$metlin_id,"</a>", sep="")), sep=""),
   #                   ifelse(match.state[i]!=1,"View",""));
    result[i, ]<-c(metabolite_names[i],
                    ifelse(match.state[i]==0, "NA", match.values[i]),
                    ifelse(match.state[i]==0, "NA", hit$hmdb_id),
                    ifelse(match.state[i]==0, "NA", hit$pubchem_id),
                    ifelse(match.state[i]==0, "NA", hit$chebi_id),
                    ifelse(match.state[i]==0, "NA", hit$kegg_id),
                    ifelse(match.state[i]==0, "NA", hit$metlin_id),
                    ifelse(match.state[i]==0, "NA", hit$smiles),
                    match.state[i]);
}
result = result[, c("Query", "Match", "HMDB", "PubChem", "ChEBI", "KEGG", "Comment")]  # 暂时仅提供这些注释内容 
write.table(result, file=output, row.names=F, sep = "\t", quote = F, na = "");

message("finish metabolite_name_map.r")

