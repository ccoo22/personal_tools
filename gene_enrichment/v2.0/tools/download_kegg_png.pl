# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Cwd qw( abs_path );
use Parallel::ForkManager;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"}; 

# 定义 -> 核心变量

 
# 检测 -> 脚本输入
my ($input_file, $output_dir, $if_help);
GetOptions(
    "input_file|i=s"  => \$input_file,
    "output_dir|o=s"  => \$output_dir,
    "help|h"          => \$if_help,
);
die help() if(defined $if_help or (not defined $input_file or not defined $output_dir));
$output_dir = Cwd::abs_path($output_dir);
###################################################################### 主程序
# 1.输出目录准备
my $html_dir = "$output_dir/html";
my $png_dir  = "$output_dir/png";
make_dir($html_dir, $png_dir);

# 2.html文件下载（文件中包含了png图路径）
download_html($input_file, $html_dir);

# 3.从html文件中，提取png下载地址
my $png_url_file = "$output_dir/png.urls.txt";
system qq{cat $html_dir/*html | grep "tmp/mark_pathway" | sed 's/<img src="//' |   sed 's/" name="pathwayimage" name="pathwayimage" usemap="#mapdata" border="0" \\/>//' |   awk '{print "http://www.genome.jp"\$0}' > $png_url_file};

# 4.下载png文件
download_png($png_url_file, $png_dir);


###################################################################### 子函数

# 下载png数据
sub download_png{
	my $png_url_file = shift @_;
	my $png_dir   = shift @_;
	
	print "下载kegg png 文件\n";
	# 提取文件中的下载地址
	my @png_addresses;
	open ADD, $png_url_file;
	while(<ADD>)
	{
		$_=~s/[\r\n]//g;
		next if($_!~/\w/);
		push @png_addresses, $_;
	}
	close ADD;

	# 下载（如果有文件下载失败，则重新循环）
	my @png_losts = @png_addresses;  # 记录还有哪些没有完成下载
	my $time = 0;
	do
	{
		$time++;
		my $max_threads = 10;
		my $pm          = Parallel::ForkManager->new($max_threads);

		# 并行下载
		foreach my $add (@png_losts) 
		{
			my $pid = $pm->start and next;
			my ($name) = $add =~ /(\w+)\.png/;
			system("wget -q -O $png_dir/$name.png '$add'  ");  # 下载
			$pm->finish;
		}
		$pm->wait_all_children;

		# 检查是否下载失败
		my @png_errors = ();
		my $lost_count = 0;
		my @lost_hsa_ids = ();
		foreach my $add (@png_losts) 
		{
			my ($name) = $add =~ /(\w+)\.png/;
			my $png_file = "$png_dir/$name.png";
			next if(-e $png_file  and is_end_png($png_file ) == 1);  # html是否完整
			push @png_errors, $add;
			push @lost_hsa_ids, $name;
			$lost_count++;
		}

		# 修改html_losts;
		@png_losts = @png_errors;
		print "[Warnings] kegg png文件下载循环 [第 $time 次]，有 [$lost_count] 个通路下载失败，重新尝试下载失败通路(睡眠10s),如果反复失败，可手工协助下载或放弃该通路 \n失败通路列表：@lost_hsa_ids\n " if($lost_count > 0);
		sleep 10;

	} while (scalar @png_losts != 0);	 
	print "kegg png 文件下载完成\n";

}


# 下载html数据
sub download_html{
	my $input_file = shift @_;
	my $html_dir   = shift @_;
	
	print "下载kegg html 文件\n";
	# 提取文件中的下载地址
	my @html_addresses;
	open ADD, $input_file;
	while(<ADD>)
	{
		$_=~s/[\r\n]//g;
		next if($_!~/\w/);
		my ($pathway_id, $add) = split /\t/, $_;
		push @html_addresses, $add;
	}
	close ADD;

	# 下载（如果有文件下载失败，则重新循环）
	my @html_losts = @html_addresses;  # 记录还有哪些没有完成下载
	my $time = 0;
	do
	{
		$time++;
		my $max_threads = 10;
		my $pm          = Parallel::ForkManager->new($max_threads);

		# 并行下载
		foreach my $add (@html_losts) 
		{
			my $pid = $pm->start and next;
			my ($name) = $add =~ /show_pathway\?(\w+)/;
			system("wget -q -O $html_dir/$name.html '$add'  ");  # 下载
			$pm->finish;
		}
		$pm->wait_all_children;

		# 检查是否下载失败
		my @html_errors = ();
		my $lost_count = 0;
		my @lost_hsa_ids = ();
		foreach my $add (@html_losts) 
		{
			my ($name) = $add =~ /show_pathway\?(\w+)/;
			my $html_file = "$html_dir/$name.html";
			next if(-e $html_file  and is_end_html($html_file ) == 1);  # html是否完整
			push @html_errors, $add;
			push @lost_hsa_ids, $name;
			$lost_count++;
		}

		# 修改html_losts;
		@html_losts = @html_errors;
		print "[Warnings] kegg html文件下载循环 [第 $time 次]，有 [$lost_count] 个通路下载失败，重新尝试下载失败通路(睡眠10s),如果反复失败，可手工协助下载或放弃该通路（下载地址长度不能超过3040个字符）\n失败通路列表：@lost_hsa_ids\n " if($lost_count > 0);
		sleep 10;

	} while (scalar @html_losts != 0);
	print "kegg html 文件下载完成\n";

}

# hmtl尾部是否有问题（存在html下载了一半，然后终止了的情况）
sub is_end_html
{
	my $html = shift;
	local $/ = undef;
	open HTML, $html;
	my $context = <HTML>;
	close HTML;
	return $context =~ /<\/html>/ ? 1 : 0;
}
# png尾部是否有问题（存在png下载了一半，然后终止了的情况）
sub is_end_png
{
    my $png = shift;
    my $buff;
    my $str = "00000000";  # 补充前缀，防止下面提取数据失败，产生报错警告
    open PNG, $png or die "Can't open $png";
    binmode(PNG);
    while(read(PNG, $buff, 1)){
        my $hex = unpack("H*", $buff);
        $str .= $hex;
    }
    close PNG; 
    my $line = substr $str, -8;

    return $line  eq qq{ae426082} ? 1 : 0;
}


# 创建目录
sub make_dir{
    my @dirs = @_;
    foreach my $dir(@dirs)
    {
        mkdir $dir if(not -e $dir);
    }
}

sub help{
    my $info = "
Program: gene_enrichment, 基因富集分析GO/KEGG
Version: 2019-12-19
Contact: 129 甘斌

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --input_file/-i   kegg通路地址文件，一行一个地址，例如: http://www.genome.jp/kegg-bin/show_pathway?hsa00052
         --output_dir/-o   结果输出路径
         --help/-h         查看帮助文档
         注意：GO分析要慢一些
    \n";
    return $info;
}
