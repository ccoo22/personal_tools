$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
基于plink家系信息的家系图绘制
Version: v1.0 2020-04-22
Contact: 129 甘斌
\n"; print $SCRIPT_INFO;
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 变量

my $DEFAULT_SOFT_MADELINE2   = "/usr/local/bin/madeline2";

# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($ped, $output, $SOFT_MADELINE2, $if_help);
GetOptions(
	"ped|p=s"           => \$ped,
	"output|o=s"        => \$output,
	"SOFT_MADELINE2=s"  => \$SOFT_MADELINE2,

	"help|h"              => \$if_help,
);
die "
Options: 必填

        --ped/-p               plink的ped文件，重点是前6列
        --output/-o            输出文件前缀

Options: 可选
        --SOFT_MADELINE2           madeline2 软件路径       (default: $DEFAULT_SOFT_MADELINE2)
        --help/-h                  查看帮助文档
\n" if (defined $if_help or not defined $ped or not defined $output);

$SOFT_MADELINE2   = $DEFAULT_SOFT_MADELINE2 if (not defined $SOFT_MADELINE2);
 
###################################################################### 主程序
my $data_madeline = "$output.madeline.data";

open PED, $ped;
open MDATA, ">$data_madeline";
my @heads = qw{FamilyId IndividualId Gender Father Mother Deceased Proband Affected};
print MDATA (join "\t", @heads) . "\n";


while(<PED>)
{
    $_=~s/[\r\n]//g;
    my ($fid, $iid, $pid, $mid, $sex, $pheno, $tmp) = split /\s+/, $_, 7;

    # 准备新的格式数据
    my ($gender, $father, $mother, $decceased, $proband, $affected) = ('.', '.', '.', '.', '.', '.');
    
    $gender = 'M' if($sex eq '1');
    $gender = 'F' if($sex eq '2');
    
    $father = $pid if($pid ne '0' and $pid ne '-9');
    $mother = $mid if($mid ne '0' and $mid ne '-9');

    $affected = 'A' if($pheno eq '2'); 
    $affected = 'U' if($pheno eq '1');

    my @values = ($fid, $iid, $gender, $father, $mother, $decceased, $proband, $affected);
    print MDATA (join "\t", @values) . "\n";
}
close PED;
close MDATA;
system("$SOFT_MADELINE2 -L 'IndividualId' $data_madeline --nolabeltruncation --color  --noiconlabels --embedded  --font OpenSans-Regular --outputprefix $output");



# open PED, $ped;
# my @heads = qw{FamilyId IndividualId Gender Father Mother Deceased Proband Affected};
# my %hashData;
# my %hashMaxLength;

# while(<PED>)
# {
#     $_=~s/[\r\n]//g;
#     my ($fid, $iid, $pid, $mid, $sex, $pheno, $tmp) = split /\s+/, $_, 7;

#     # 准备新的格式数据
#     my ($gender, $father, $mother, $decceased, $proband, $affected) = ('.', '.', '.', '.', '.', '.');
    
#     $gender = 'M' if($sex eq '1');
#     $gender = 'F' if($sex eq '2');
    
#     $father = $pid if($pid ne '0' and $pid ne '-9');
#     $mother = $mid if($mid ne '0' and $mid ne '-9');

#     $affected = 'A' if($pheno eq '2');
#     $affected = 'U' if($pheno eq '1');

#     $hashData{$.}{'FamilyId'}     = $fid;
#     $hashData{$.}{'IndividualId'} = $iid;
#     $hashData{$.}{'Gender'}       = $gender;
#     $hashData{$.}{'Father'}       = $father;
#     $hashData{$.}{'Mother'}       = $mother;
#     $hashData{$.}{'Deceased'}     = $decceased;
#     $hashData{$.}{'Proband'}      = $proband;
#     $hashData{$.}{'Affected'}     = $affected;


#     # 记录每一列最长的字符长度
#     foreach my $title(keys %{$hashData{$.}})
#     {
#         my $length = length($hashData{$.}{$title});
#         $hashMaxLength{$title} = $length if(not exists $hashMaxLength{$title} or $hashMaxLength{$title} < $length);
#     }
# }
# close PED;

# # 数据输出
# open MDATA, ">$data_madeline";

# # 表头部分
# map{ print MDATA "$_\n" } @heads;
# print MDATA "\n";

# # 数据部分
# foreach my $row(sort {$a<=>$b} keys %hashData)
# {   
#     my $string = '';
#     foreach my $title(@heads)
#     {   
#         my $value      = $hashData{$row}{$title};
#         my $length     = length($value);
#         my $max_length = $hashMaxLength{$title};
#         my $extend     = $max_length - $length + 1;  # 需要补充空格的数量
#         my $extend_blank = ' ' x $extend;

#         $string .= "$value$extend_blank";
#     }
#     print MDATA "$string\n";
# }


# close MDATA;

# system("$SOFT_MADELINE2 -L 'IndividualId' $data_madeline --nolabeltruncation --color  --noiconlabels --embedded  --outputprefix $output");

###################################################################### 子程序


