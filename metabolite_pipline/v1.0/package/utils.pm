package utils;

use strict;
use warnings;
use Exporter;
use vars qw(@EXPORT_OK @ISA);

@ISA = 'Exporter';
@EXPORT_OK = qw(is_file_ok is_pdf_ok check_result_file parse_group);

# 检验文件是否为空
sub is_file_ok{
    my @files = @_;
    my $isOK = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file == 0);
    }    
    return $isOK;
}
# 检验pdf文件是否为空
sub is_pdf_ok{
    my @files = @_;
    my $isOK = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file < 4400);
    }    
    return $isOK;
}

sub check_result_file{
    my @results          = @_;

    my $file_type = "";
    my $is_ok     = 1;  # 默认全部正常
    my @errors;
    foreach my $result(@results)
    {   
        # 文件类型确认
        if($result eq 'txt' or $result eq 'pdf')
        {
            $file_type = $result;
            next;
        }

        # 检查txt文件
        if($file_type eq 'txt' and is_file_ok($result) == 0)
        {
            $is_ok = 0;
            push @errors, $result;
        }

        # 检查pdf文件
        if($file_type eq 'pdf' and is_pdf_ok($result) == 0)  # pdf文件现在先用file模式检测，后面增加pdf专用检测方案
        {
            $is_ok = 0;
            push @errors, $result;
        }
    }

    my $error = join ", ", @errors;
    return $error;
}


sub parse_group{
    my $group_names   = shift @_;
    my $sample_counts = shift @_;
    my $group_num = scalar(@$group_names);  # 分组总数量

    my %hashGroup;

    # 对比构成
	foreach my $index(0..$#$group_names)
    {   
        # （1）所有分组一起
        my $group_name_1 = $$group_names[$index];
        $hashGroup{'all_samples'}{'GroupCount'}++;  # 分组数量
        $hashGroup{'all_samples'}{'SampleCount'} += $$sample_counts[$index];  # 样本数量

        # （2）两两对比
        if($group_num > 2 and $index != $#$group_names)
        {
            foreach my $index2(($index + 1)..$#$group_names)
            {   
                my $group_name_2 = $$group_names[$index2];
                $hashGroup{"$group_name_1\_VS_$group_name_2"}{'GroupCount'}  = 2;  # 分组数量
                $hashGroup{"$group_name_1\_VS_$group_name_2"}{'SampleCount'} = $$sample_counts[$index] + $$sample_counts[$index2];  # 样本数量
                # $hashGroup{"$group_name_1\_VS_$group_name_2"}{'SampleCount1'} = $$sample_counts[$index];  # 样本数量
                # $hashGroup{"$group_name_1\_VS_$group_name_2"}{'SampleCount2'} = $$sample_counts[$index2];  # 样本数量
            }
        }
	}
    return %hashGroup;
}



1