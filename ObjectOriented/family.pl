
use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_registry_from_url('mysql://anonymous@ensembldb.ensembl.org:5306/');

my $member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Member');
my $member = $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE','ENSG00000105697');

my $family_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','Family');
my $families = $family_adaptor->fetch_all_by_Member($member);

foreach my $family (@{$families}) {
#    print join " ", map { $family->$_ }  qw(description description_score),"\n";
    print $family->description . " " . $family->description_score . "\n";


    foreach my $member_attribute (@{$family->get_all_Member_Attribute}) {
        my ($member, $attribute) = @{$member_attribute};
        print $member->stable_id," ",$member->taxon_id,"\n";
    }

    my $simple_align = $family->get_SimpleAlign();
    my $alignIO = Bio::AlignIO->newFh(
        -interleaved => 0,
        -file => ">simple_align.phy",
        -format => "phylip",
        -idlength => 20);

    print $alignIO $simple_align;

    $simple_align = $family->get_SimpleAlign('cdna');
    my $alignIO = Bio::AlignIO->newFh(
        -interleaved => 0,
        -file => ">cdna_align.phy",
        -format => "phylip",
        -idlength => 20);

    print $alignIO $simple_align;
}
