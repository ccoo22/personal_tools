# 导入 -> 系统 package
use strict;
use warnings;
use Encode;
use File::Spec;
use Getopt::Long;
use Excel::Writer::XLSX;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 检测 -> 脚本输入
my ($input_file, $output_file, $sheet_name, $if_help);
GetOptions(
    "sheet_name|n=s" => \$sheet_name,
    "input_file|i=s" => \$input_file,
    "output_file|o=s" => \$output_file,
    "help|h" => \$if_help,
);
die help() if(defined $if_help or (not defined $input_file or not defined $output_file or not defined $sheet_name));
###################################################################### 主程序
my $workbook = Excel::Writer::XLSX->new($output_file);
my %format   = table_format($workbook);

my @file_lists = split /;/, $input_file;
my @name_lists = split /;/, $sheet_name;

# 如果表格名称大于31字符，则需要替换
my @sheet_name_maps;

foreach my $count(0..$#file_lists)
{
    my $file       = $file_lists[$count];
    my ($name, @marks) = split /,/, $name_lists[$count];

    # 需要标记P值的注释信息
    my $mark_count = scalar(@marks) / 2;
    my %hashPvalue = map{ my $tmp = 2 * $_; ($marks[$tmp], $marks[$tmp + 1])   } 0..$mark_count-1;
    print "process $file \n";

    # 创建表格对象
    my $length_limit = 31;
    if(length($name) > $length_limit)
    {   # 名字太长，替换
        my $prefix = "Sheet_$count\_";
        my $name_new  = $prefix . substr($name, 0, $length_limit - length($prefix));
        my @maps = ($name_new, $name);
        push @sheet_name_maps, \@maps;

        $name = $name_new;
    }
    my $sheet = $workbook->add_worksheet($name);
    $sheet->set_row(0, 60);
    
    open INPUT, $file;

    # 表头处理
    my $line1 = <INPUT>;
       $line1 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line1;
    
    my $row = 0;
    my @formats = map{ 'title'} 0..$#heads;
    output_sheet(\@heads, \@formats, $row, $sheet, \%format, $name);
    $row++;

    # 正文处理
    my $is_excess_line = 0; # 行数是否超过了excel表的限制
    while(<INPUT>)
    {
        $_=~s/[\r\n]//g;
        next if($_!~/\w/);
        my @tmps = split /\t/, $_;

        # 数据处理，定义颜色格式
        my @datas;
        my @formats;
        foreach my $col(0..$#heads)
        {   
            my $head  = $heads[$col];
            my $value = (exists $tmps[$col]) ? $tmps[$col] : ""; 
            
            my $format_define = 'normal'; # 表格默认格式normal
               $format_define = 'orange' if(exists $hashPvalue{$head} and $hashPvalue{$head} =~ /\d/ and $value =~ /\d/ and $value < $hashPvalue{$head}); # Pvalue < defined

            push @datas, $value;
            push @formats, $format_define;
        }
        output_sheet(\@datas, \@formats, $row, $sheet, \%format, $name);
        $row++;
        
        # excel表有行数限制
        if($row >= 1048576)
        {   
            $is_excess_line = 1;
            last;
        }
    }
    close INPUT;

    if($is_excess_line ==1)
    {
        print "[Warnings] file line is excess 1048576, excel can not hold so many lines. so we only output 1048576 line to excel 【$file】\n";
    }
}

# 输出名字映射关系
if(scalar(@sheet_name_maps) > 0)
{   
    print "[Warnings] 因为表格名称过31个字符，脚本自动重命名表格名称，并生成映射关系表格\n";
    my $sheet = $workbook->add_worksheet('Sheet_Map');
    my $row = 0;
    foreach my $map_index(@sheet_name_maps)
    {   
        my @formats = map{ 'normal' }  0..$#$map_index;
        output_sheet($map_index, \@formats, $row, $sheet, \%format, 'Sheet_Map');
        $row++;
    }
}



###################################################################### 子函数

sub output_sheet{
    my $datas   = shift @_; # 数据
    my $formats = shift @_; # 数据对应的格式
    my $row     = shift @_; # 行
    my $sheet   = shift @_; # 表格对象
    my $format  = shift @_; # 格式数据库列表
    my $name    = shift @_; # 表格名称，主要是应对readme的中文问题

    foreach my $col(0..$#$datas)
    {   
        $$datas[$col] = decode("gb2312", $$datas[$col]) if($name =~ /readme/i); # 中文解码
        $sheet->write($row, $col, $$datas[$col], $format->{$$formats[$col]});
    }
}

sub table_format{
    my ($workbook)=@_;
    my %format=();
    $format{'title'} = $workbook->add_format();
    $format{'title'} ->set_align('center');
    $format{'title'} ->set_align('vcenter');
    $format{'title'} ->set_size(12);
    $format{'title'} ->set_font("Times New Roman");
    $format{'title'} ->set_border();
    $format{'title'} ->set_bg_color("yellow");
    $format{'title'} ->set_color("black");
    
    $format{'normal'} = $workbook->add_format();
    $format{'normal'} ->set_align('center');
    $format{'normal'} ->set_align('vcenter');
    $format{'normal'} ->set_size(12);
    $format{'normal'} ->set_font("Times New Roman");
    $format{'normal'} ->set_border();
    
    $format{'small'} = $workbook->add_format();
    $format{'small'} ->set_align('vcenter');
    $format{'small'} ->set_size(10);
    $format{'small'} ->set_font("Times New Roman");
    $format{'small'} ->set_border();
    
    $format{'seq'} = $workbook->add_format();
    $format{'seq'} ->set_align('vcenter');
    $format{'seq'} ->set_size(11);
    $format{'seq'} ->set_font("Courier New");
    $format{'seq'} ->set_border();
    
    $format{'left'} = $workbook->add_format();
    $format{'left'} ->set_align('vcenter');
    $format{'left'} ->set_size(12);
    $format{'left'} ->set_font("Times New Roman");
    $format{'left'} ->set_border();
    
    $format{'orange'} = $workbook->add_format();
    $format{'orange'} ->set_align('vcenter');
    $format{'orange'} ->set_size(12);
    $format{'orange'} ->set_font("Times New Roman");
    $format{'orange'} ->set_bg_color("#fac090");
    $format{'orange'} ->set_border();

    $format{'skyblue'} = $workbook->add_format();
    $format{'skyblue'} ->set_align('vcenter');
    $format{'skyblue'} ->set_size(12);
    $format{'skyblue'} ->set_font("Times New Roman");
    $format{'skyblue'} ->set_bg_color("#538ed5");
    $format{'skyblue'} ->set_border();

    $format{'bold'} = $workbook->add_format( bold => 1 );
    $format{'blue'} = $workbook->add_format( color => "#538ed5" );
    $format{'redbold'} = $workbook->add_format( color => "#ff0000", bold => 1, );
    $format{'italic'} = $workbook->add_format( italic => 1 );
    $format{'boldblue'} = $workbook->add_format( bold => 1, color => "#538ed5" );
    $format{'bolditalic'} = $workbook->add_format( bold => 1, italic => 1 );
    $format{'blueitalic'} = $workbook->add_format( color => "#538ed5", italic => 1 );
    $format{'boldblueitalic'} = $workbook->add_format( bold => 1, color => "#538ed5", italic => 1 );
    
    $format{'readme5'} = $workbook->add_format();
    $format{'readme5'}->set_align('vcenter');
    $format{'readme5'}->set_size(11);
    $format{'readme5'}->set_font("Times New Roman");
    $format{'readme5'}->set_border();
    $format{'readme5'}->set_text_wrap();

    
    $format{'readme1'} = $workbook->add_format();
    $format{'readme1'}->set_align('center');
    $format{'readme1'}->set_align('vcenter');
    $format{'readme1'}->set_bold();
    $format{'readme1'}->set_size(14);
    $format{'readme1'}->set_font("Times New Roman");
    $format{'readme1'}->set_border();

    $format{'readme2'} = $workbook->add_format();
    $format{'readme2'}->set_align('vcenter');
    $format{'readme2'}->set_bold();
    $format{'readme2'}->set_size(14);
    $format{'readme2'}->set_font("Times New Roman");

    $format{'readme2tmp'} = $workbook->add_format();
    $format{'readme2tmp'}->set_right();

    $format{'readme3'} = $workbook->add_format();
    $format{'readme3'}->set_align('center');
    $format{'readme3'}->set_align('vcenter');
    $format{'readme3'}->set_bold();
    $format{'readme3'}->set_size(11);
    $format{'readme3'}->set_font("Times New Roman");
    $format{'readme3'}->set_border();

    $format{'readme4'} = $workbook->add_format();
    $format{'readme4'}->set_align('vcenter');
    $format{'readme4'}->set_bold();
    $format{'readme4'}->set_size(11);
    $format{'readme4'}->set_font("Times New Roman");
    $format{'readme4'}->set_border();

    $format{'readme5'} = $workbook->add_format();
    $format{'readme5'}->set_align('vcenter');
    $format{'readme5'}->set_size(11);
    $format{'readme5'}->set_font("Times New Roman");
    $format{'readme5'}->set_border();
    $format{'readme5'}->set_text_wrap();

    return %format;
}


sub help{
    my $info = "
Program: convert table to excel ##### 'tab' seperate table #####
Version: 2018-12-31
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --input_file/-i    輸入文件，例如text.txt;text2.txt 可支持多个文件，“分号”分割,sheet名字与-s一一对应
         --output_file/-o   輸出文件, 例如restult.xlsx
         --sheet_name/-s    表格名稱, 例如：FPKM,ABC,0.05,CDE,0.01;FPKM2,AAA,0.01;FPKM3
                                       注：支持指定列P值按阈值标颜色,可以同时写多列，每两列为一个单位，用“逗号”分割
         --help/-h          查看帮助文档
    \n";
    return $info;
}

