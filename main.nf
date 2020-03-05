#!/usr/bin/env nextflow

/*
 * Copyright (c) 2017-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'XXXXXX'.
 *
 *   XXXXXX is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   XXXXXX is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with XXXXXX.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main XXX pipeline script
 *
 * @authors
 * Edgar Garriga

 */

/*
 * defaults parameter definitions
 */


//singularity run docker://evolbioinfo/pastml --tree ~/CBCRG/nf_pastML/test/tree.nwk --data ~/CBCRG/nf_pastML/test/encoded.csv --columns c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22 c23 c24 c25 c26 c27 c28 c29 c30 c31 c32 c33 c34 c35 c36 c37 c38 c39 c40 c41 c42 c43 c44 c45 c46 c47 c48 c49 c50 c51 c52 c53 --html_compressed ~/CBCRG/nf_pastML/test/Albanian_map.html --data_sep ,

// input sequences to align in fasta format
params.alignments = "$baseDir/test/*.aln"
//params.metadata = "$baseDir/test/*.csv"
params.trees = "$baseDir/test/*.dnd"
params.separator= ","
params.predictionMethod="JOINT"

// output directory
//defined in nextflow.config

log.info """\
         F  A  M  S  A    P  i  p  e  l  i  n  e  ~  version 0.1"
         ======================================="
         Input alignments (FASTA)                       : ${params.alignments}
         Input trees (NEWICK)                           : ${params.trees}
         Prediction Method                              : ${params.predictionMethod}
         Output directory (DIRECTORY)                   : ${params.outdir}
         """
         .stripIndent()

// Channels containing sequences
if ( params.alignments ) {
  Channel
  .fromPath(params.alignments)
  .map { item -> [ item.baseName.tokenize('.')[0], item] }
  .view()
  .set { aln }
}
/**
if ( params.metadata ) {
  Channel
  .fromPath(params.metadata)
  .map { item -> [ item.baseName, item] }
  .set { meta }
}**/

// Channels for user provided trees
if ( params.trees ) {
  Channel
    .fromPath(params.trees)
    .map { item -> [ item.baseName.tokenize('.')[0], item] }
    .view()
    .set { trees }
}

process createMetadata{
  conda 'environment.yml'
  tag "metadata-${id}"
  publishDir "${params.outdir}/metadata", mode: 'copy', overwrite: true

  input:
    set val(id), file(alignment) from aln

  output:
    set val(id), file("${id}.csv"), file("command.txt") into metaOut

  script:
  """
  echo "\$(pastml.py ${alignment} ${id})"
  """

}

trees
  .combine ( metaOut, by:0 )
  .into { treesAndMeta; treesAndMeta2 }

treesAndMeta2.view()

process pastML {
    tag "pastML-${id}"
    publishDir "${params.outdir}/pastML", mode: 'copy', overwrite: true

    input:
    set val(id), file(tree), file(metadata), file(command) from treesAndMeta

    output:
    file("${id}_out.html") into pastmlOut

    script:
    """
    var="\$(cat ${command})"
    pastml --tree ${tree} --data ${metadata} --data_sep ${params.separator} --columns \$var --prediction_method ${params.predictionMethod} --html_compressed ${id}_out.html 
    """
}

workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}
