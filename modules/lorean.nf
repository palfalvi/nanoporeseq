process lorean {

  process.container = 'docker://lfaino/lorean'
  //process.containerOptions '-B <PATH_TO_AUGUSTUS_CONF_FOLDER>:/opt/LoReAn/third_party/software/augustus/config/ -B <PATH_TO_LIBRARY_FOLDER>:/usr/local/RepeatMasker/Libraries/'
  singularity.enabled = true

  publishDir "${params.outdir}/lorean", mode: 'copy'

  input:
    path genome
    path protein_ref

  output:
    path "*", emit: annotation

  script:

    def prefix      = params.lorean_prefix   ? "--prefix_gene ${params.prefix}"       : ""
    def stranded    = params.lorean_stranded ? "--stranded"                           : ""
    def iproscan    = params.lorean_iproscan ? "--interproscan"                       : ""
    def adapters    = params.lorean_adapters ? "--adapter ${params.lorean_adapters}"  : ""
    def long_reads  = params.lorean_long     ? "--long_reads ${params.lorean_long}"   : ""
    def short_reads = params.lorean_short    ? "--short_reads ${params.lorean_short}" : ""
    def species     = params.lorean_species  ? "--species ${params.lorean_species}"   : "--species Xx"

    """
    sleep 10
    lorean --help > lorean.test
    """
//    """
//    wget https://github.com/lfaino/LoReAn/raw/master/third_party/software/config.augustus.tar.gz && tar -zxvf config.augustus.tar.gz

//    wget https://github.com/lfaino/LoReAn/raw/master/third_party/software/RepeatMasker.Libraries.tar.gz && tar -xvzf RepeatMasker.Libraries.tar.gz

//    lorean \
//    --threads $task.cpus \
//    --minimap2 \
//    --max_intron_length 10000 \
//    -pr $protein_ref \
//    -sp $species \
//    $long_reads \
//    $short_reads \
//    $prefix \
//    $stranded \
//    $iproscan \
//    $adapter \
//    $genome
//    """
}
