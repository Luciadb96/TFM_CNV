import os, sys
import gzip
import configparser, optparse

######################################################################

class OptionParser(optparse.OptionParser):

	def check_required(self, opt):
		option = self.get_option(opt)
		atrib = getattr(self.values, option.dest)
		
		if atrib is None:
			return False
		else:
			return True

###########################################################################
		
def Set_Chr_Nr_ (Chr):
    """ Sort by chromosome """
    if Chr:
        New = Chr[3:]
        if New == 'X': New = 23
        elif New == 'Y': New = 24
        elif New == 'M': New = 25
        else: New = int(New)
    else:
        New = 0
    return New
	
###########################################################################

def run(argv=None):
	
	if argv is None: argv = sys.argv  ##lista de argumentos de la línea de comandos pasados a un script de Python. argv[0] es el nombre del script

	parser = OptionParser(add_help_option=True,description="The script parses variant call files of estd/nstd studies from dbVar")
	parser.add_option("--i",default=None,help="tsv.gz input file",dest="input_filename")
	parser.add_option("--o",default=None,help="output path",dest="output_path")
	
	(options, args) = parser.parse_args(argv[1:])

	if len(argv) == 1:
		sys.exit(0)

	if not parser.check_required("--i"):
		raise IOError('The input file has not been provided')
	if not parser.check_required("--o"):
		raise IOError('The output path has not been provided')
	
	filename_tsv = options.input_filename
	if not os.path.exists(filename_tsv):
		raise IOError('The input file provided does not exist. %s' % (filename_tsv))
	
	output_path = options.output_path
	if not os.path.exists(output_path):
		raise IOError('The output path provided does not exist. %s' % (output_path))

	fi = gzip.open(filename_tsv, 'rt') ##with gzip.open(filename_tsv,'rb') as fi:
	l_lines = map(lambda x: x.strip().split('\t'), fi.readlines())
	fi.close()
	
	hash_cnvs = {}
	
	"""
		Fields expected in tsv.gz file:
	
		0.variant_call_accession
		1.variant_call_id
		2.variant_call_type
		3.experiment_id
		4.sample_id
		5.sampleset_id
		6.assembly
		7.chr
		8.contig
		9.outer_start
		10.start
		11.inner_start
		12.inner_stop
		13.stop
		14.outer_stop
		15.insertion_length
		16.variant_region_acc
		17.variant_region_id
		18.copy_number
		19.description
		20.validation
		21.zygosity
		22.origin
		23.phenotype
		24.hgvs_name
		25.placement_method
		26.placement_rank
		27.placements_per_assembly
		28.remap_alignment
		29.remap_best_within_cluster
		30.remap_coverage
		31.remap_diff_chr
		32.remap_failure_code
		
		Specifications found in: https://github.com/ncbi/dbvar/blob/master/specs/dbVar.xsd
		
		Filters are applyied in order to get a consistent and non-redundant set of variant according to https://github.com/ncbi/dbvar/blob/master/Structural_Variant_Sets/Nonredundant_Structural_Variants/README.md
		
		Filters_
		
		- 'origin': Only variant calls are from germline samples only (no somatic)
		- placements are "BestAvailable" on the assembly (guarantees no duplicate placements for a variant)
		- placements are on finished chromosomes only (not on NT_ or NW_ contigs)
		- placements are 1-based in the .tsv files
		- placements are zero-based start and 1-based stop in .bed and .bedpe files
		- insertions submitted to dbVar without insertion_length or submitted sequence are not included in the NR files
		
	"""

	for line in l_lines:

		if line[0][0] == '#':
			continue

		rank = line[26]
		if rank != 'BestAvailable':  # This is the best placement available in the assembly
			continue

		origin = line[22]
		if origin == "somatic":
			continue

		contig = line[8]
		if contig.find("NT_") != -1 or contig.find("NW_") != -1:  ##find() method returns the index of first occurrence of the substring (if found). if not found, it returns -1
			continue

		chrom = line[7]
		chrom = "chr%s" % (chrom)

		sample = line[4]
		variant_id = line[0]
		variant_region_acc = line[16]
		variant_region_id = line[17]
		cnv_type = line[2].replace(" ", "_")
		start = int(line[11])
		end = int(line[12])
		number_cp = line[18]
		if number_cp == "":
			number_cp = '.'
			
		# zygosity = line[21]

		num_map = line[27]
		if num_map == "":
			num_map = '.'

		validation = line[20]
		if validation == "":
			validation = '.'

		hash_cnvs.setdefault(sample, []).append(
			(chrom, min(start, end), max(start, end), cnv_type, number_cp, variant_id, variant_region_acc, variant_region_id, validation, num_map))
		
	l_samples = sorted(hash_cnvs.keys())
	
	for i in l_samples:
		
		l_cnv = hash_cnvs[i] 
		
		l_cnv_sorted = sorted(l_cnv, key=lambda x: (Set_Chr_Nr_(x[0]), x[1], x[2]))
		
		fo = open(os.path.join(output_path, "%s_CNV_all.bed" % (i)), 'w')
		fo.write("#chrom\tstart\tend\tcnv_type\tnumber_cp\tvariant_id\tvariant_region_acc\tvariant_region_id\tvalidation\tnum_map\n")
		fo.write("\n".join(map(lambda x: "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9]), l_cnv_sorted)))
		fo.close()


###################################################################################################################	

if __name__=='__main__':
	try:
		run()
	except Exception as e: 
		print(e)
