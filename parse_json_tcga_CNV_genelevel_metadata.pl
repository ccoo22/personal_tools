# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use JSON;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量


# 检测 -> 脚本输入
my ($input_json, $output, $if_help);
GetOptions(
    "input_json|i=s"   => \$input_json,
    "output|o=s"       => \$output,
    "help|h"           => \$if_help,
);
die help() if(defined $if_help or (not defined $input_json or not defined $output));
###################################################################### 主程序

my $json_string = read_json_string($input_json);
open OUTPUT, ">$output";
print OUTPUT "sample_name\tcase_id\tentity_id\n";

my $json = new JSON;
my $json_obj = $json->decode($json_string);
foreach my $obj_one(@$json_obj)
{
    foreach my $obj_associated_entities(@{$obj_one->{'associated_entities'}})
    {
        my $sample_name  = $obj_associated_entities->{'entity_submitter_id'};
        my $case_id      = $obj_associated_entities->{'case_id'};
        my $entity_id    = $obj_associated_entities->{'entity_id'};
        print OUTPUT "$sample_name\t$case_id\t$entity_id\n";  
    }
}


close OUTPUT;



###################################################################### 子函数

sub read_json_string{
    my $fson_file = shift @_;

    my $json_string = '';
    open FILTERJSON ,$fson_file;
    while(<FILTERJSON>)
    {
        $json_string.=$_;
    }
    close FILTERJSON;
    return $json_string;
}


sub help{
    my $info = "
Program: parse JSON metadata of TCGA CNV GeneLevel ， 拆分TCGA CNV GeneLevel下载数据json格式的metadata文件
Version: 2019-03-11
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --input_json/-i    json文件
         --output/-o        输出文件
         --help/-h          查看帮助文档
    \n";
    return $info;
}

