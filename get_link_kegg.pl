use strict;
use warnings;
use File::Spec;
use Getopt::Long;
$|=1;

 
my ($kegg, $count, $output_dir, $preffix, $help);
GetOptions(
	"kegg|k=s"        => \$kegg,
	"count|c=s"       => \$count,
	"output_dir|o=s"  => \$output_dir,
	"preffix|p=s"     => \$preffix,
	"help|h!"	      => \$help,
);
die help() if (defined $help or (not defined $kegg or not defined $output_dir or not defined $preffix));
$count = 10 if(not defined $count);
############################################ 主流程
mkdir $output_dir if(not -e $output_dir);
my %hashKEGG = read_kegg($kegg, $count);

my $link = "$output_dir/$preffix.link.txt";
my $node = "$output_dir/$preffix.node.txt";
open NODE, ">$node";
open LINK, ">$link";
print LINK "source\ttarget\n";
print NODE "node\tvalue\tcolor\tsize\n";
foreach my $id(sort keys %hashKEGG)
{   
    my $desc = $hashKEGG{$id}{'Description'};
    print NODE "$desc\thsa\tblue\t30\n";
    my @genes = split /\//, $hashKEGG{$id}{'geneID'};
    foreach my $gene(@genes)
    {
        print LINK "$desc\t$gene\n";
        print NODE "$gene\tgene\tred\t30\n";
    }
}
close LINK;
close NODE;

##############子程序##############

sub read_kegg{
    my $kegg = shift @_;
    my $count = shift @_;

    my %hashKEGG;
    open FILE, $kegg;
    my $line1 = <FILE>;
       $line1 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line1;

    my $ok = 0;
    while(<FILE>)
    {
        $_ =~ s/[\r\n]//g;
        next if($_!~/\w/);
        my @datas = split /\t/, $_;
        my %hashTmp;
        foreach my $col(0..$#heads)
        {
            my $value = exists $datas[$col] ? $datas[$col] : "";
            $hashTmp{$heads[$col]} = $value;
        }

        my $pvalue = $hashTmp{'pvalue'};
        my $id     = $hashTmp{'ID'};
        next if($pvalue >= 0.05);
        
        map{ $hashKEGG{$id}{$_} = $hashTmp{$_} } keys %hashTmp;
         $ok++;
        last if($ok >= $count);  # 只要前count个
    }
    close FILE;
    return %hashKEGG;


}

sub help{
    my $info = "
Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
        --kegg/-k          kegg分析表格
        --count/-c         选取前n个做图。默认： 10
        --output_dir/-o    输出路径
        --preffix/-p       输出文件前缀
        --help/-h          查看帮助文档
    \n";
    return $info;  
}

