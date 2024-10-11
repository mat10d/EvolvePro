import os
from evolvepro.process.exp_process import *

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))

# Define the WT sequences for all proteins
wt_sequences = {
    't7_pol': 'MNTINIAKNDFSDIELAAIPFNTLADHYGERLAREQLALEHESYEMGEARFRKMFERQLKAGEVADNAAAKPLITTLLPKMIARINDWFEEVKAKRGKRPTAFQFLQEIKPEAVAYITIKTTLACLTSADNTTVQAVASAIGRAIEDEARFGRIRDLEAKHFKKNVEEQLNKRVGHVYKKAFMQVVEADMLSKGLLGGEAWSSWHKEDSIHVGVRCIEMLIESTGMVSLHRQNAGVVGQDSETIELAPEYAEAIATRAGALAGISPMFQPCVVPPKPWTGITGGGYWANGRRPLALVRTHSKKALMRYEDVYMPEVYKAINIAQNTAWKINKKVLAVANVITKWKHCPVEDIPAIEREELPMKPEDIDMNPEALTAWKRAAAAVYRKDKARKSRRISLEFMLEQANKFANHKAIWFPYNMDWRGRVYAVSMFNPQGNDMTKGLLTLAKGKPIGKEGYYWLKIHGANCAGVDKVPFPERIKFIEENHENIMACAKSPLENTWWAEQDSPFCFLAFCFEYAGVQHHGLSYNCSLPLAFDGSCSGIQHFSAMLRDEVGGRAVNLLPSETVQDIYGIVAKKVNEILQADAINGTDNEVVTVTDENTGEISEKVKLGTKALAGQWLAYGVTRSVTKRSVMTLAYGSKEFGFRQQVLEDTIQPAIDSGKGLMFTQPNQAAGYMAKLIWESVSVTVVAAVEAMNWLKSAAKLLAAEVKDKKTGEILRKRCAVHWVTPDGFPVWQEYKKPIQTRLNLMFLGQFRLQPTINTNKDSEIDAHKQESGIAPNFVHSQDGSHLRKTVVWAHEKYGIESFALIHDSFGTIPADAANLFKAVRETMVDTYESCDVLADFYDQFADQLHESQLDKMPALPAKGNLNLRDILESDFAFA',
    'r2': 'VKVTVPDKNPPCPCCSTRLNSVLALIDHLKGSHGKRRVCFRCAKCGRENFNHHSTVCHFAKCKGPSEEKPPVGEWICEVCGRDFTTKIGLGQHKRLAHPMVRNQERIDASQPKETSNRGAHKKCWTKEEEELLARLEVQFEGHKNINKLIAEHITTKTNKQISDKRRQMTRKDKGEGGAAGKLGPDTGRGNHSQAKVGNNGLGGNQLPGGPAATKDKAGCHLDKEEGNRIAISQQKKGRLQGRYHKEIKRRLEEGVINTFTKAFKQLLECQEVQPLINKTAQDCFGLLESACHIRTALRGKNKKETQEKPTGGQCLKWMKKRAVKKGNYLRFQRLFHLDRGKLARIILDDIECLSCDIAPSEIYSVFKARWETPGQFAGLGNFKSTGKADNKAFSDLITAKEIKKNVQEMSKGSAPGPDGIAIGDIKGMDPGYSRTAELFNLWLTSGEIPDMVRGCRTVLIPKSTQPERLKDINNWRPITIGSILLRLFSRIITARMTKACPLNPRQRGFIRAAGCSENLKLLQTIIRTAKSEHRPLGVVFVDIAKAFDTVSHQHILHVLQQRGVDPHIIGLVSNMYKDISTFVTTKKDTHTDKIQIRVGVKQGDPLSPLLFNLAMDPLLCKLEESGNGFHRGGHTITAMAFADDLVLLSDSWENMEKNIEILEAFCDLTGLKTQGQKCHGFYIKPTKDSYTVNNCAAWTIYGTPLNMINPGDSEKYLGLQIDPWTGIARSNISSKLDSWLERINQAPLKPLQKLDILKTYTIPRLTYMVDHSEMKAGALEALDLQIRSAVKDWLHLPSCTCDAILYVSTKDGGLGVTKLAGLIPSIQARRLHRIAQSPDETMKAFLDKEQMEKQYAKLWVQAGGKREKIPSIWDALPTPVLLTTSDTLSEWEAPNPKSKYPRPCNWRRKEFEKWTKLQCQGRGIQNFKGDVISNNWIQNYRRIPHRKLLTAVQLRANVYPTREFLGRGRGDDCVKFCRHCEVDLETCGHIISYCPVTKEARIKRHNRICERLIEEAEKKDWVVFKEPHIRDAVKELFKPDLIFVKEDRALVVDVTVRFEATTTSLEEAAIEKVDKYKRLETEVRSLTNAKDVLFMGFPLGARGKWYQGNFKLLDMLGLSESRQVTVAKTLSTDALISSVDIVHMFASKARKMNLVTV',
    'fanzor': 'MKRKREDLTLWDAANVHKHKSMWYWWEYIRRKDMVNHEKTDCDVIQLLQSASVKKQKTQSDKFLTSFSVGIRPTKHQKRVLNEMLRVSNYTYNWCLWLVNEKGLKPHQFELQKIVCKTNANDVDPQYRMENDDWFFNNKMTSVKLTSCKNFCTSYKSAKSLKSKLKRPMSVSNIIQGSFCVPKLFIRHLSSKDVSTDNTNMQNRYICMMPDNFEKRSNPKERFLKLAKPITKIPPIDHDVKIVKRADGMFIMNIPCDPKYTRRNASNDTIEKRVCGIDPGGRTFATVYDPIDCCVFQVGIKEDKQYVISKLHNKIDHAHMHLTKAQNKKQQQAARERIVSLKKTHLKLKTFVDDIHLKLSSHLVKEYQYVALGKINVAQLVKTDRPKPLSKRAKRDLLYWQHYRFRQRLTHRTTNTECILDVQNEAYTSKTCGVCGTINKNLEKSETFYCDQCKYNTHRDVNGARNILLKSLRMFPFEKQQQ',
    'mlv': 'TLNIEDEYRLHETSKEPDVSLGSTWLSDFPQAWAETGGMGLAVRQAPLIIPLKATSTPVSIKQYPMSQEARLGIKPHIQRLLDQGILVPCQSPWNTPLLPVKKPGTNDYRPVQDLREVNKRVEDIHPTVPNPYNLLSGLPPSHQWYTVLDLKDAFFCLRLHPTSQPLFAFEWRDPEMGISGQLTWTRLPQGFKNSPTLFNEALHRDLADFRIQHPDLILLQYVDDLLLAATSELDCQQGTRALLQTLGNLGYRASAKKAQICQKQVKYLGYLLKEGQRWLTEARKETVMGQPTPKTPRQLREFLGKAGFCRLFIPGFAEMAAPLYPLTKPGTLFNWGPDQQKAYQEIKQALLTAPALGLPDLTKPFELFVDEKQGYAKGVLTQKLGPWRRPVAYLSKKLDPVAAGWPPCLRMVAAIAVLTKDAGKLTMGQPLVILAPHAVEALVKQPPDRWLSNARMTHYQALLLDTDRVQFGPVVALNPATLLPLPEEGLQHNCLDILAEAHGTRPDLTDQPLPDADHTWYTDGSSLLQEGQRKAGAAVTTETEVIWAKALPAGTSAQRAELIALTQALKMAEGKKLNVYTDSRYAFATAHIHGEIYRRRGWLTSEGKEIKNKDEILALLKALFLPKRLSIIHCPGHQKGHSAEARGNRMADQAARKAAITETPDTSTLLIENSSP',
    'bxb1': 'MRALVVIRLSRVTDATTSPERQLESCQQLCAQRGWDVVGVAEDLDVSGAVDPFDRKRRPNLARWLAFEEQPFDVIVAYRVDRLTRSIRHLQQLVHWAEDHKKLVVSATEAHFDTTTPFAAVVIALMGTVAQMELEAIKERNRSAAHFNIRAGKYRGSLPPWGYLPTRVDGEWRLVPDPVQRERILEVYHRVVDNHEPLHLVAHDLNRRGVLSPKDYFAQLQGREPQGREWSATALKRSMISEAMLGYATLNGKTVRDDDGAPLVRAEPILTREQLEALRAELVKTSRAKPAVSTPSLLLRVLFCAVCGEPAYKFAGGGRKHPRYRCRSMGFPKHCGNGTVAMAEWDAFCEEQVLDLLGDAERLEKVWVAGSDSAVELAEVNAELVDLTSLIGSPAYRAGSPQREALDARIAALAARQEELEGLEARPSGWEWRETGQRFGDWWREQDTAAKNTWLRSMNVRLTFDVRGGLTRTIDFGDLQEYEQHLRLGSVVERLHTGMS',
    'ca': 'MTVTDDYLANNVDYASGFKGPLPMPPSKHIAIVACMDARLDVYRMLGIKEGEAHVIRNAGCVVTDDVIRSLAISQRLLGTREIILLHHTDCGMLTFTDDDFKRAIQDETGIRPTWSPESYPDAVEDVRQSLRRIEVNPFVTKHTSLRGFVFDVATGKLNEVTPAAALEARKEAELAAATAEQ',
    'mmfunc': 'MKRKREQMTLWKAAFVNGQETFKSWIDKARMLELNCDVSSASSTHYSDLNLKTKCAKMEDKFMCTFSVGIRPTSKQKRTLNQMLKVSNHAYNWCNYLVKEKDFKPKQFDLQRVVTKTNSTDVPAEYRLPGDDWFFDNKMSSIKLTACKNFCTMYKSAQTNQKKTKVDLRNKDIAMLREGSFEVQKKYVRLLTEKDIPDERIRQSRIALMADNFSKSKKDWKERFLRLSKNVSKIPPLSHDMKVCKRPNGKFVLQIPCDPIYTRQIQVHTSDSICSIDPGGRTFATCYDPSNIKAFQIGPEADKKEVIHKYHEKIDYVHRLLAYAQKKKQTQAVQDRIGQLKKLHLKLKTYVDDVHLKLCSYLVKNYKLVVLGKISVSSIVRKDRPNHLAKKANRDLLCWQHYRFRQRLLHRVRGTDCEAIAQDERYTSKTCGNCGVKNNKLGGKETFICESCNYKTHRDVNGARNILCKYLGLFPFAA',
    'psacas12f': 'MPSETYITKTLSLKLIPSDEEKQALENYFITFQRAVNFAIDRIVDIRSSFRYLNKNEQFPAVCDCCGKKEKIMYVNISNKTFKFKPSRNQKDRYTKDIYTIKPNAHICKTCYSGVAGNMFIRKQMYPNDKEGWKVSRSYNIKVNAPGLTGTEYAMAIRKAISILRSFEKRRRNAERRIIEYEKSKKEYLELIDDVEKGKTNKIVVLEKEGHQRVKRYKHKNWPEKWQGISLNKAKSKVKDIEKRIKKLKEWKHPTLNRPYVELHKNNVRIVGYETVELKLGNKMYTIHFASISNLRKPFRKQKKKSIEYLKHLLTLALKRNLETYPSIIKRGKNFFLQYPVRVTVKVPKLTKNFKAFGIDRGVNRLAVGCIISKDGKLTNKNIFFFHGKEAWAKENRYKKIRDRLYAMAKKLRGDKTKKIRLYHEIRKKFRHKVKYFRRNYLHNISKQIVEIAKENTPTVIVLEDLRYLRERTYRGKGRSKKAKKTNYKLNTFTYRMLIDMIKYKAEEAGVPVMIIDPRNTSRKCSKCGYVDENNRKQASFKCLKCGYSLNADLNAAVNIAKAFYECPTFRWEEKLHAYVCSEPDK'
}

# Generate single amino acid mutants for t7_pol
t7_pol_wt_fasta = os.path.join(project_root, 'data', 'exp', 'wt_fasta', 't7_pol_WT.fasta')
generate_wt(wt_sequences['t7_pol'], t7_pol_wt_fasta)
t7_pol_output_file = os.path.join(project_root, 'output', 'exp', 't7_pol.fasta')
generate_single_aa_mutants(t7_pol_wt_fasta, t7_pol_output_file)

# Generate single amino acid mutants for r2
r2_wt_fasta = os.path.join(project_root, 'data', 'exp', 'wt_fasta', 'r2_WT.fasta')
generate_wt(wt_sequences['r2'], r2_wt_fasta)
r2_output_file = os.path.join(project_root, 'output', 'exp', 'r2.fasta')
generate_single_aa_mutants(r2_wt_fasta, r2_output_file)

# Generate single amino acid mutants for fanzor
fanzor_wt_fasta = os.path.join(project_root, 'data', 'exp', 'wt_fasta', 'fanzor_WT.fasta')
generate_wt(wt_sequences['fanzor'], fanzor_wt_fasta)
fanzor_output_file = os.path.join(project_root, 'output', 'exp', 'fanzor.fasta')
generate_single_aa_mutants(fanzor_wt_fasta, fanzor_output_file)

# Generate single amino acid mutants for mlv
mlv_wt_fasta = os.path.join(project_root, 'data', 'exp', 'wt_fasta', 'mlv_WT.fasta')
generate_wt(wt_sequences['mlv'], mlv_wt_fasta)
mlv_output_file = os.path.join(project_root, 'output', 'exp', 'mlv.fasta')
generate_single_aa_mutants(mlv_wt_fasta, mlv_output_file)

# Generate single amino acid mutants for bxb1
bxb1_wt_fasta = os.path.join(project_root, 'data', 'exp', 'wt_fasta', 'bxb1_WT.fasta')
generate_wt(wt_sequences['bxb1'], bxb1_wt_fasta)
bxb1_output_file = os.path.join(project_root, 'output', 'exp', 'bxb1.fasta')
generate_single_aa_mutants(bxb1_wt_fasta, bxb1_output_file)

# Generate single amino acid mutants for ca
ca_wt_fasta = os.path.join(project_root, 'data', 'exp', 'wt_fasta', 'ca_WT.fasta')
generate_wt(wt_sequences['ca'], ca_wt_fasta)
ca_output_file = os.path.join(project_root, 'output', 'exp', 'ca.fasta')
generate_single_aa_mutants(ca_wt_fasta, ca_output_file)

# Generate single amino acid mutants for mmfunc
mmfunc_wt_fasta = os.path.join(project_root, 'data', 'exp', 'wt_fasta', 'mmfunc_WT.fasta')
generate_wt(wt_sequences['mmfunc'], mmfunc_wt_fasta)
mmfunc_output_file = os.path.join(project_root, 'output', 'exp', 'mmfunc.fasta')
generate_single_aa_mutants(mmfunc_wt_fasta, mmfunc_output_file)

# Generate single amino acid mutants for psacas12f
psacas12f_wt_fasta = os.path.join(project_root, 'data', 'exp', 'wt_fasta', 'psacas12f_WT.fasta')
generate_wt(wt_sequences['psacas12f'], psacas12f_wt_fasta)
psacas12f_output_file = os.path.join(project_root, 'output', 'exp', 'psacas12f.fasta')
generate_single_aa_mutants(psacas12f_wt_fasta, psacas12f_output_file)

# Generate n-mutant combinations for t7_pol for 2 to 7 mutations
t7_pol_wt_fasta = os.path.join(project_root, 'data', 'exp', 'wt_fasta', 't7_pol_WT.fasta')
t7_pol_mutant_file = os.path.join(project_root, 'data', 'exp', 'n_mutant_dicts', 't7_pol_n_mutants.xlsx')
t7_pol_2_output_file = os.path.join(project_root, 'output', 'exp', 't7_pol_2nd.fasta')
t7_pol_3_output_file = os.path.join(project_root, 'output', 'exp', 't7_pol_3rd.fasta')
t7_pol_4_output_file = os.path.join(project_root, 'output', 'exp', 't7_pol_4th.fasta')
t7_pol_5_output_file = os.path.join(project_root, 'output', 'exp', 't7_pol_5th.fasta')
t7_pol_6_output_file = os.path.join(project_root, 'output', 'exp', 't7_pol_6th.fasta')
t7_pol_7_output_file = os.path.join(project_root, 'output', 'exp', 't7_pol_7th.fasta')

generate_n_mutant_combinations(t7_pol_wt_fasta, t7_pol_mutant_file, n=2, output_file=t7_pol_2_output_file, threshold=1)
generate_n_mutant_combinations(t7_pol_wt_fasta, t7_pol_mutant_file, n=3, output_file=t7_pol_3_output_file, threshold=1)
generate_n_mutant_combinations(t7_pol_wt_fasta, t7_pol_mutant_file, n=4, output_file=t7_pol_4_output_file, threshold=1)
generate_n_mutant_combinations(t7_pol_wt_fasta, t7_pol_mutant_file, n=5, output_file=t7_pol_5_output_file, threshold=1)
generate_n_mutant_combinations(t7_pol_wt_fasta, t7_pol_mutant_file, n=6, output_file=t7_pol_6_output_file, threshold=1)
generate_n_mutant_combinations(t7_pol_wt_fasta, t7_pol_mutant_file, n=7, output_file=t7_pol_7_output_file, threshold=1)

# Generate n-mutant combinations for mlv for 2 to 7 mutations
mlv_wt_fasta = os.path.join(project_root, 'data', 'exp', 'wt_fasta', 'mlv_WT.fasta')
mlv_mutant_file = os.path.join(project_root, 'data', 'exp', 'n_mutant_dicts', 'mlv_n_mutants.xlsx')
mlv_2_output_file = os.path.join(project_root, 'output', 'exp', 'mlv_2nd.fasta')
mlv_3_output_file = os.path.join(project_root, 'output', 'exp', 'mlv_3rd.fasta')
mlv_4_output_file = os.path.join(project_root, 'output', 'exp', 'mlv_4th.fasta')
mlv_5_output_file = os.path.join(project_root, 'output', 'exp', 'mlv_5th.fasta')
mlv_6_output_file = os.path.join(project_root, 'output', 'exp', 'mlv_6th.fasta')
mlv_7_output_file = os.path.join(project_root, 'output', 'exp', 'mlv_7th.fasta')

generate_n_mutant_combinations(mlv_wt_fasta, mlv_mutant_file, n=2, output_file=mlv_2_output_file, threshold=1)
generate_n_mutant_combinations(mlv_wt_fasta, mlv_mutant_file, n=3, output_file=mlv_3_output_file, threshold=1)
generate_n_mutant_combinations(mlv_wt_fasta, mlv_mutant_file, n=4, output_file=mlv_4_output_file, threshold=1)
generate_n_mutant_combinations(mlv_wt_fasta, mlv_mutant_file, n=5, output_file=mlv_5_output_file, threshold=1)
generate_n_mutant_combinations(mlv_wt_fasta, mlv_mutant_file, n=6, output_file=mlv_6_output_file, threshold=1)
generate_n_mutant_combinations(mlv_wt_fasta, mlv_mutant_file, n=7, output_file=mlv_7_output_file, threshold=1)

# Generate n-mutant combinations for mlv_2 for 2 to 7 mutations
mlv_wt_fasta = os.path.join(project_root, 'data', 'exp', 'wt_fasta', 'mlv_WT.fasta')
mlv_mutant_file = os.path.join(project_root, 'data', 'exp', 'n_mutant_dicts', 'mlv_n_mutants_2.xlsx')
mlv_2_output_file = os.path.join(project_root, 'output', 'exp', 'mlv_2nd_2.fasta')
mlv_3_output_file = os.path.join(project_root, 'output', 'exp', 'mlv_3rd_2.fasta')
mlv_4_output_file = os.path.join(project_root, 'output', 'exp', 'mlv_4th_2.fasta')
mlv_5_output_file = os.path.join(project_root, 'output', 'exp', 'mlv_5th_2.fasta')
mlv_6_output_file = os.path.join(project_root, 'output', 'exp', 'mlv_6th_2.fasta')
mlv_7_output_file = os.path.join(project_root, 'output', 'exp', 'mlv_7th_2.fasta')

generate_n_mutant_combinations(mlv_wt_fasta, mlv_mutant_file, n=2, output_file=mlv_2_output_file, threshold=1)
generate_n_mutant_combinations(mlv_wt_fasta, mlv_mutant_file, n=3, output_file=mlv_3_output_file, threshold=1)
generate_n_mutant_combinations(mlv_wt_fasta, mlv_mutant_file, n=4, output_file=mlv_4_output_file, threshold=1)
generate_n_mutant_combinations(mlv_wt_fasta, mlv_mutant_file, n=5, output_file=mlv_5_output_file, threshold=1)
generate_n_mutant_combinations(mlv_wt_fasta, mlv_mutant_file, n=6, output_file=mlv_6_output_file, threshold=1)
generate_n_mutant_combinations(mlv_wt_fasta, mlv_mutant_file, n=7, output_file=mlv_7_output_file, threshold=1)