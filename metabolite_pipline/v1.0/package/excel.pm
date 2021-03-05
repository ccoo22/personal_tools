package excel;

use strict;
use warnings;
use Exporter;
use vars qw(@EXPORT_OK @ISA);

@ISA = 'Exporter';
@EXPORT_OK = qw(set_excel_format txt_to_excel);

sub set_excel_format{
	my $workbook = shift;
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

	$format{'bold'} = $workbook->add_format( bold => 1 );
	$format{'blue'} = $workbook->add_format( color => "#538ed5" );
	$format{'redbold'} = $workbook->add_format( color => "#ff0000", bold => 1, );
	$format{'italic'} = $workbook->add_format( italic => 1 );
	$format{'boldblue'} = $workbook->add_format( bold => 1, color => "#538ed5" );
	$format{'bolditalic'} = $workbook->add_format( bold => 1, italic => 1 );
	$format{'blueitalic'} = $workbook->add_format( color => "#538ed5", italic => 1 );
	$format{'boldblueitalic'} = $workbook->add_format( bold => 1, color => "#538ed5", italic => 1 );

	$format{'readme_title'} = $workbook->add_format();
	$format{'readme_title'}->set_align('vcenter');
	$format{'readme_title'}->set_bold();
	$format{'readme_title'}->set_size(24);
	$format{'readme_title'}->set_font("Times New Roman");

	$format{'readme_bold'} = $workbook->add_format();
	$format{'readme_bold'}->set_align('vcenter');
	$format{'readme_bold'}->set_bold();
	$format{'readme_bold'}->set_size(12);
	$format{'readme_bold'}->set_font("Times New Roman");
	$format{'readme_bold'}->set_border();
	$format{'readme_bold'}->set_text_wrap();

	$format{'readme'} = $workbook->add_format();
	$format{'readme'}->set_align('vcenter');
	$format{'readme'}->set_size(12);
	$format{'readme'}->set_font("Times New Roman");
	$format{'readme'}->set_border();
	$format{'readme'}->set_text_wrap();

	return %format;
}


sub txt_to_excel{
    my $workbook = shift @_;
    my $format   = shift @_;
    my $file     = shift @_;
    my $name     = shift @_;

    my $sheet = $workbook->add_worksheet($name);
    open FILE, $file;
    my $line1 = <FILE>;
       $line1 =~s/[\r\n]//g;
    my @heads = split /\t/, $line1;

    my $row = 0;
    foreach my $col(0..$#heads)
    {
        $sheet->write($row, $col, $heads[$col], $format->{"title"});
    }

    $row++;
    while(<FILE>)
    {
        $_=~s/[\r\n]//g;
        my @datas = split /\t/, $_;
        foreach my $col(0..$#heads)
        {   
            my $value = exists $datas[$col] ? $datas[$col] : "";
            $sheet->write($row, $col, $value, $format->{"normal"});
        }
        $row++;
    }
    close FILE;
    
}


1