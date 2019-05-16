

/*


 --Calibrated ChIP-Seq Workflow--
      -for paired-end reads-
	   version 1.0


 Usage - 'nextflow run calibrated.nf <parameter specifications>'


 // software versions \\


 - FastQC 0.11.2
 - MultiQC 1.7
 - BWA 0.7.17
 - Samtools 1.3.1
 - genomeCoverageBed 2.27.0
 - wigToBigWig 377-0
 - MACSII 2.1.2

 
 // pipeline input parameters \\


 - params.reads_all takes all read files in defined directory. All used read files should be included under this parameter.
 - params.reads_IP is the immunoprecipitated reads
 - params.reads_input is the input reads
 - params.forward_primer is the forward primer sequence
 - params.reverse_primer is the reverse primer sequence
 - params.genome_calib is the fasta sequence of the calibration organism genome
 - params.genome_exper is the fasta sequence of the experimental organism genome
 - params.chrom_size is the chromosome sizes text file of the experimental genome
	
*/

params.reads_all ="$baseDir/fastq/input/APR_18_VM/*.fastq.gz"
params.reads_IP = "$baseDir/fastq/input/APR_18_VM/IP-449*R{1,2}*.fastq.gz"
params.reads_input = "$baseDir/fastq/input/APR_18_VM/Input-449*R{1,2}*.fastq.gz"
params.forward_primer ='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
params.reverse_primer ='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
params.genome_calib = "$baseDir/genomes/W303.fa"
params.genome_exper = "$baseDir/genomes/pombe.fa"
params.chrom_size = "$baseDir/genomes/pombe.chrom.sizes"

/* channel creation of parameters */

Channel.fromPath(params.reads_all).set{ reads_file }		

genome_calib_file = file(params.genome_calib)
genome_exper_file = file(params.genome_exper)


forward = file(params.forward_primer)
reverse = file(params.reverse_primer)


Channel.fromFilePairs(params.reads_IP, flat: true) 
       .set{read_pairs_IP}					


Channel.fromFilePairs(params.reads_input, flat: true)
       .set{read_pairs_input}

chrom_size = file(params.chrom_size)  
     


process fastqc {
	publishDir 'fastqc'
	tag "FastQC on $reads_file_in"

	input:
	file reads_file_in from reads_file

	output:
	file '*_fastqc.{zip,html}' into fastqc_o
		

	script:
	"""
	
	fastqc -q "$reads_file_in"   
	"""
}

process multiqc {
	tag "multiqc"
	publishDir "fastqc"
	
	input:
	file (multiqc_in:'multiqc_in/*') from fastqc_o.collect()

	output:
	file '*'

	script:
	"""
	multiqc *
	cp multiqc_report.html /homes/matts74/public_html/project/multiqc
	
	"""
	
}


process index {
	tag "Indexing on $genome_calib_file and $genome_exper_file"
	publishDir 'genome_index'	

	input:
	
	file genome_calib from genome_calib_file
	file genome_exper from genome_exper_file

	output:
	file "BWAindex" into genome_index_ch		

	script:
	"""
	bwa index "$genome_calib"
	mkdir BWAindex && mv ${genome_calib}* BWAindex 
	bwa index "$genome_exper"
	mv ${genome_exper}* BWAindex
	"""
}



process trim_nf {
	tag "Trim on $read_pairs"
	publishDir 'trim_nf'
		
	input:
	set val(id), file(read1), file(read2) from read_pairs
	val forward from forward
	val reverse from reverse

	output:
        file "*.fq.gz" into trimmed_reads
	file "*trimming_report.txt" into trimgalore_results
	file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports	


	script:
	"""	
	mkdir "$id"

	trim_galore --paired -a '$forward' -a2 '$reverse' -q 20 --length 65 --fastqc '$read1' '$read2'
	
	""" 
}



process bwa_IP_first_run {
	
	tag "IP alignment on $id to experimental and calibration"
	

	input:

	set val(id), file(bwa_IP_first_run_read_1), file(bwa_IP_first_run_read_2) from read_pairs_IP	
	file index from genome_index_ch	

	output:
	file '*exper.bam' into bwa_IP_first_run_exper_bam
	file '*calib.bam' into bwa_IP_first_run_calib_bam
	file '*exper.bam' into macs_IP_bam_in	

	script:
	"""
	bwa mem ${index}/pombe.fa '$bwa_IP_first_run_read_1' '$bwa_IP_first_run_read_2'| samtools view -b > '${id}_exper.bam'
	bwa mem ${index}/W303.fa '$bwa_IP_first_run_read_1' '$bwa_IP_first_run_read_2'| samtools view -b > '${id}_calib.bam'

	"""
}

process bwa_IP_second_run_calib {
	
	tag "IP unmapped alignment to calibration"
        publishDir 'aligned_bwa_bam'

        input:


	file bwa_IP_second_run_calib_bam_in from bwa_IP_first_run_exper_bam
        file index from genome_index_ch

        output:
        file '*IPci.bam' into IPci_bam

        script:
        """
        samtools view -bh -f 4 '$bwa_IP_second_run_calib_bam_in' > bwa_IP_second_run_calib_bam_unmapped.bam
	samtools sort -n -@ 5 bwa_IP_second_run_calib_bam_unmapped.bam -o bwa_IP_second_run_calib_unmapped_sort.bam
	samtools bam2fq bwa_IP_second_run_calib_unmapped_sort.bam > bwa_IP_second_run_calib_unmapped_sort.fastq
	
	bwa mem ${index}/W303.fa bwa_IP_second_run_calib_unmapped_sort.fastq | samtools view -b > ${bwa_IP_second_run_calib_bam_in.baseName}_IPci.bam

        """
}


process bwa_IP_second_run_exper {

        tag "IP unmapped alignment to experiment"
        publishDir 'aligned_bwa_bam'

        input:

        file bwa_IP_second_run_exper_bam_in from bwa_IP_first_run_calib_bam
        file index from genome_index_ch

        output:
        file '*IPxi.bam' into IPxi_bam

        script:
        """
        samtools view -bh -f 4 '$bwa_IP_second_run_exper_bam_in' > bwa_IP_second_run_exper_bam_unmapped.bam
        samtools sort -n -@ 5 bwa_IP_second_run_exper_bam_unmapped.bam -o bwa_IP_second_run_exper_unmapped_sort.bam
        samtools bam2fq bwa_IP_second_run_exper_unmapped_sort.bam > bwa_IP_second_run_exper_unmapped_sort.fastq

        bwa mem ${index}/pombe.fa bwa_IP_second_run_exper_unmapped_sort.fastq | samtools view -b > ${bwa_IP_second_run_exper_bam_in.baseName}_IPxi.bam

        """
}

process bwa_input_first_run {

        tag "input alignment to experimental and calibration"


        input:

        set val(id), file(bwa_input_first_run_read_1), file(bwa_input_first_run_read_2) from read_pairs_input
        file index from genome_index_ch

        output:
        file '*exper.bam' into bwa_input_first_run_exper_bam
        file '*calib.bam' into bwa_input_first_run_calib_bam
	file '*exper.bam' into macs_input_bam_in
        script:
        """
        bwa mem ${index}/pombe.fa '$bwa_input_first_run_read_1' '$bwa_input_first_run_read_2'| samtools view -b > '${id}_exper.bam'
        bwa mem ${index}/W303.fa '$bwa_input_first_run_read_1' '$bwa_input_first_run_read_2'| samtools view -b > '${id}_calib.bam'

        """
}

process bwa_input_second_run_calib {

        tag "input unmapped alignment to calibration"
        publishDir 'aligned_bwa_bam'

        input:

        file bwa_input_second_run_calib_bam_in from bwa_input_first_run_exper_bam
        file index from genome_index_ch

        output:
        file '*Wci.bam' into Wci_bam

        script:
        """
        samtools view -bh -f 4 '$bwa_input_second_run_calib_bam_in' > bwa_input_second_run_calib_bam_unmapped.bam
        samtools sort -n -@ 5 bwa_input_second_run_calib_bam_unmapped.bam -o bwa_input_second_run_calib_unmapped_sort.bam
        samtools bam2fq bwa_input_second_run_calib_unmapped_sort.bam > bwa_input_second_run_calib_unmapped_sort.fastq

        bwa mem ${index}/W303.fa bwa_input_second_run_calib_unmapped_sort.fastq | samtools view -b > ${bwa_input_second_run_calib_bam_in.baseName}_Wci.bam

        """
}

process bwa_input_second_run_exper {

        tag "input unmapped alignment to experiment"
        publishDir 'aligned_bwa_bam'

        input:

        file bwa_input_second_run_exper_bam_in from bwa_input_first_run_calib_bam
        file index from genome_index_ch

        output:
        file '*Wxi.bam' into Wxi_bam

        script:
        """
        samtools view -bh -f 4 '$bwa_input_second_run_exper_bam_in' > bwa_input_second_run_exper_bam_unmapped.bam
        samtools sort -n -@ 5 bwa_input_second_run_exper_bam_unmapped.bam -o bwa_input_second_run_exper_unmapped_sort.bam
        samtools bam2fq bwa_input_second_run_exper_unmapped_sort.bam > bwa_input_second_run_exper_unmapped_sort.fastq

        bwa mem ${index}/pombe.fa bwa_input_second_run_exper_unmapped_sort.fastq | samtools view -b > ${bwa_input_second_run_exper_bam_in.baseName}_Wxi.bam

        """
}


process samtools_processing {
	
	tag "post-alignment processing"
	publishDir 'samtools_o'

	input:
	
	file samtool_bam_in from IPxi_bam

	output:
	
	set file('*.sorted.bam'), file('*sorted.bam.bai') into ch_samStats
	file '*.stats.txt' into stats_out

	script:
	"""
	samtools sort -@ 5 '$samtool_bam_in' -o ${samtool_bam_in.baseName}.sorted.bam
	samtools index ${samtool_bam_in.baseName}.sorted.bam ${samtool_bam_in.baseName}.sorted.bam.bai
	samtools stats ${samtool_bam_in.baseName}.sorted.bam > ${samtool_bam_in.baseName}.stats.txt  

	"""

} 

process RPM_num {

	input:
	
	file RPM_no from stats_out

	output:
        
	stdout into reads_no_out
	stdout into calib_reads_no_out

        script:

        """
        cat '$RPM_no' | grep 'reads mapped:' | cut -f 3 

        """

}

process RPM {
		
	input:
	val RPM_no from reads_no_out.toInteger()

	output:
	stdout into RPM_out

	shell:
	'''
	echo 1000000/$((!{RPM_no})) | bc -l
	'''
}

process uncalib_bedgraph {

	tag "bedgraph"
	publishDir 'bedgraph_o'


	input:
	set file(sorted_bam_in), file(bam_index) from ch_samStats
	val num_value from RPM_out

	output:
	file '*_rpm.bedgraph' into bedgraph_out

	script:
	"""
	genomeCoverageBed -ibam -bg -i '$sorted_bam_in' -scale '$num_value' -bg  > ${sorted_bam_in.baseName}_uncalib_rpm.bedgraph
        cp *.bedgraph /homes/matts74/public_html/project/bedgraph
	"""

}

bedgraph_out.into { bedgraph_out1; bedgraph_out2 }

/*
process uncalib_bigwig {

	tag "bigwig"
	publishDir 'bigwig'

	input:
	file bedgraph from bedgraph_out1
	file chrom_size_no from chrom_size 
	
	output:
	file '*.bw' into uncalib_bigwig_out

	script:
	"""
	wigToBigWig '$bedgraph' '$chrom_size_no' ${bedgraph.baseName}_uncalib.bw
	"""

}

*/
process Wxi_calculation {
	
	
	input:

	file Wxi_bam from Wxi_bam
	
	output:
	stdout Wxi_num into Wxi_no

	shell:
	'''	

	samtools stats !{Wxi_bam} > Wxi_no.stats
	cat Wxi_no.stats | grep "reads mapped:" | cut -f 3
	'''
}


process Wci_calculation {


        input:

        file Wci_bam from Wci_bam

        output:
        stdout Wci_num into Wci_no

        shell:
        '''
	samtools stats !{Wci_bam} > Wci_no.stats
        cat Wci_no.stats | grep "reads mapped:" | cut -f 3
	'''
}


process IPci_calculation {


        input:

        file IPci_bam from IPci_bam

        output:
        stdout IPci_num into IPci_no

        shell:
	'''
	samtools stats !{IPci_bam} > IPci_no.stats
        cat IPci_no.stats | grep "reads mapped:" | cut -f 3
	'''
}

process WciIPxi {
	
	input:
	val Wci_no from Wci_no
	val IPxi_no from calib_reads_no_out

	output:
	stdout into WciIPxi_out

	script:
	"""
	perl -e "print $Wci_no*$IPxi_no" 
	"""
}


process WxiIPci {

        input:
        val Wxi_no from Wxi_no
        val IPci_no from IPci_no

        output:
        stdout into WxiIPci_out

        script:
        """
        perl -e  "print $Wxi_no*$IPci_no"
	"""
}


	
process inver_WxiIPci_calculation {

	input:
	val WxiIPci from WxiIPci_out

	output:
	stdout into inver_WxiIPci

	script:
	"""
	echo '1/$WxiIPci'| bc -l 
	
	"""
}

process OR_calculation { 
	
	input:
	val inver_WxiIPci from inver_WxiIPci
	val WciIPxi from WciIPxi_out

	output:
	stdout into OR_val
	
	script:
	"""
	perl -e "print $WciIPxi*$inver_WxiIPci" 
	"""
}


process calibrated_bedgraph {
	
	"calibration of bedgraph"
	publishDir 'bedgraph_o'

	input:
	file calib_in from bedgraph_out2
	val OR from OR_val

	output:
	file '*.bedgraph' into calib_bedgraph

	shell:
	'''
	less !{calib_in} | awk '{print $1 "\t" $2 "\t" $3 "\t" $4*!{OR}}' > !{calib_in.baseName}_calibrated.bedgraph 
	cp *calibrated.bedgraph ~/public_html/project
	'''
}
/*
process calibrated_bigWig {
	
	"calibrated bigWig"
	publishDir 'bigwig'
	
	input:
	file calib_bedgraph from calib_bedgraph

	output:
	file '*calib.bw'
	
	script:
	"""
	wigToBigWig '$calib_bedgraph' '$chrom_size_no_2' ${calib_bedgraph.baseName}_calib.bw        
	"""
}

*/
process MACSII {

	tag "Peak calling"
	publishDir 'macs_o'

	input:
	file IP_bam from macs_IP_bam_in
	file input_bam from macs_input_bam_in

	output:
	file '*.{bed,r,xls,bdg}' into macs_out
	file '*.narrowPeak' into peaks_out

	script:
	"""
	macs2 callpeak -B --nomodel -n '${IP_bam.baseName}' -c '$input_bam' -t '$IP_bam' -q 0.01 -g 1.2e7
	"""
}

