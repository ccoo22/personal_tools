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
print OUTPUT "sample_name\tcase_id\tfile_id\tfile_name\n";

my $json = new JSON;
my $json_obj = $json->decode($json_string);
foreach my $obj_one(@$json_obj)
{
    my $file_name  = $obj_one->{'file_name'};
    my $file_id  = $obj_one->{'file_id'};
    my $obj_entity = $obj_one->{'associated_entities'};
    foreach my $obh_entity_one(@$obj_entity)
    {
        my $sample_name = $obh_entity_one->{'entity_submitter_id'};
        my $case_id     = $obh_entity_one->{'case_id'};
        print OUTPUT "$sample_name\t$case_id\t$file_id\t$file_name\n";
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
Program: parse JSON metadata of TCGA FPKM ， 拆分TCGA FPKM下载数据json格式的metadata文件
         该文件可用于后续识别样本名，合并数据
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

