with open(config['SAMPLES']) as fp:
    samples = fp.read().splitlines()

rule all:
         input:
            #Prepare samples
            #================
            expand("galore/{sample}_R1_001_val_1.fq.gz", sample = samples),
            expand("galore/{sample}_R2_001_val_2.fq.gz", sample = samples),
            expand("{sample}.sam", sample = samples),
            expand("{sample}.bam", sample = samples), 
            #expand("{sample}.bam", sample =samples),
            #expand("{sample}.sorted.bam", sample =samples),
            #expand("{sample}.sorted.rmDup.bam", sample =samples),
            #expand("{sample}.bigwig", sample = samples),
            
            #expand("macs/{sample}_summits.bed", sample = samples),
            #expand("macs/{sample}_peaks.narrowPeak", sample = samples),
            #expand("Motif_{sample}/seq.autonorm.tsv", sample = samples),           
            #expand("{sample}.annotatednarrowpeaks", sample = samples),
            #expand("{sample}.annotatednarrowpeaks.stats", sample = samples),
            #expand("{sample}.annotatednarrowpeaks.bed", sample = samples),
            #expand("{sample}.bed", sample = samples),
            #expand("{sample}.annotatednarrowpeaks.fasta", sample = samples)

rule trim: 
       input: 
           r1 = "{sample}_R1_001.fastq.gz",
           r2 = "{sample}_R2_001.fastq.gz"
       output: 
           "galore/{sample}_R1_001_val_1.fq.gz",
           "galore/{sample}_R2_001_val_2.fq.gz"
       shell: 
           """
           mkdir -p galore 
           mkdir -p fastqc 
           trim_galore --gzip --retain_unpaired --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
           """ 
rule align:
              input:               
                   "galore/{sample}_R1_001_val_1.fq.gz",
                   "galore/{sample}_R2_001_val_2.fq.gz"
              params:
                   index=config['INDEX'],
                   mem = config['MEMORY'],
                   cores = config['CORES']
              output:
                   "{sample}.sam",
                   "{sample}_hist.txt" 
              shell:
                   """
                   bowtie2 --local --very-sensitive-local --no-mixed --no-discordant --phred33 -I 10 -X 700 -p {params.cores} -x {params.index} -1 {input[0]} -2 {input[1]} -S {output[0]}  &> {output[1]}
                   """
rule samTobam:
             input: 
                 "{sample}.sam",
             output: 
                 "{sample}.bam"
             shell: 
                   """
                   samtools view -bS {input} > {output}
                   """


rule merge: 
            input: 
                 "{sample}_L001.bam",
                 "{sample}_L002.bam",
                 "{sample}_L003.bam",
                 "{sample}_L004.bam"
            output:
                    "{sample}.bam"
            shell:
                   """
                   samtools merge {output} {input[0]} {input[1]} {input[2]} {input[3]} 
                   """
rule sort: 
             input: 
                  "{sample}.bam"
             output: 
                    "{sample}.sorted.bam"
             shell: 
                   """ 
                   picard SortSam I={input}  O={output} SORT_ORDER=coordinate 
                   """ 

rule remove_duplicates:
       input: 
        "{sample}.sorted.bam"
       output: 
         "{sample}.sorted.rmDup.bam",
         "{sample}.rmDup.txt"
       shell: 
           """
            picard MarkDuplicates I={input} O={output[0]} REMOVE_DUPLICATES=true METRICS_FILE={output[1]} 
           """

rule index: 
      input: 
         "{sample}.sorted.rmDup.bam"
      output: 
         "{sample}.sorted.rmDup.bam.bai"
      shell: 
          """
          samtools index {input} 
          """ 
rule bamCoverage: 
       input: 
        "{sample}.sorted.rmDup.bam",
        "{sample}.sorted.rmDup.bam.bai" 
       output: 
        "{sample}.bigwig" 
       params: 
         genome_size = config['Genome_Size'], 
         binsize = config['BINSIZE'], 
         num_processors = config['Num_Processors'] 
       shell: 
          """ 
          bamCoverage -b {input[0]} -p {params.num_processors}  --normalizeUsing RPGC --effectiveGenomeSize {params.genome_size} --binSize {params.binsize} -o {output} 
          """ 

rule macs_bed: 
      input: 
         "{sample}.sorted.rmDup.bam"
      params: 
         "{sample}", 
         genome_size = config['Genome_Size'] 
      output: 
          "macs/{sample}_summits.bed",
          "macs/{sample}_peaks.narrowPeak", 
      shell: 
           """
           bash macs2.sh 
           """


rule annotateNarrowPeaks: 
      input: 
         "macs/{sample}_peaks.narrowPeak" 
      params: 
           genome= config['GENOME'], 
           gtf = config['GTF']  
      output: 
         "{sample}.annotatednarrowpeaks", 
         "{sample}.annotatednarrowpeaks.stats"
      shell: 
          """
          annotatePeaks.pl {input} {params.genome} -gtf {params.gtf}   -annStats {output[1]}  > {output[0]}    
          """ 


rule findMotifs:
      input:
          "macs/{sample}_summits.bed"
      params:
          genome = config['GENOME'],
          output_dir = "Motif_{sample}"
      output:
          "Motif_{sample}/seq.autonorm.tsv"
      shell:
         """
          findMotifsGenome.pl {input} {params.genome} {params.output_dir} -size 200 -mask
         """


rule toBed: 
    input: 
        "{sample}.annotatednarrowpeaks"
    output: 
        "{sample}.annotatednarrowpeaks.bed"
    shell: 
        "cut -f 1-6 {input} > {output}" 


rule adjustBed: 
     input: 
         "{sample}.annotatednarrowpeaks.bed"
     output: 
         "{sample}.bed" 
     shell: 
        "sed '1d' {input} | cut -f2- > {output} " 
     
rule getFasta: 
    input: 
       "{sample}.bed"
    params: 
        genome = config['GENOME'], 
    output:
        "{sample}.annotatednarrowpeaks.fasta" 
    shell: 
        """
         bedtools getfasta -fi {params.genome} -bed {input} -fo {output}
        """
