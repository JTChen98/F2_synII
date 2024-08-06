## prepare transcripts.fa
/opt/working/projects/prj_054_1000seq_RNAseq/gffread/gffread /opt/working/projects/prj_054_1000seq_RNAseq/data/BY4741.gff \
														  -g /opt/working/projects/prj_054_1000seq_RNAseq/data/BY4741.fsa \
														  -w /opt/working/projects/prj_054_1000seq_RNAseq/data/BY4741.transcripts.fa

/opt/working/projects/prj_054_1000seq_RNAseq/gffread/gffread /opt/working/projects/prj_054_1000seq_RNAseq/data/synII_new.gff \
														  -g /opt/working/projects/prj_054_1000seq_RNAseq/data/synII_new.fsa \
														  -w /opt/working/projects/prj_054_1000seq_RNAseq/data/synII_new.transcripts.fa

## creat Xmatrix
/opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/bin/TxIndexer -t /opt/working/projects/prj_054_1000seq_RNAseq/data/BY4741.transcripts.fa \
																			 -o /opt/working/projects/prj_054_1000seq_RNAseq/TxIndexer_idx_BY

Rscript /opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/R/genPolyesterSimulation.R /opt/working/projects/prj_054_1000seq_RNAseq/data/BY4741.transcripts.fa /opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_BY

/opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/bin/GenTC -i /opt/working/projects/prj_054_1000seq_RNAseq/TxIndexer_idx_BY \
																		 -l IU -1 /opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_BY/sample_01_1.fasta \
																		 -2 /opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_BY/sample_01_2.fasta \
																		 -o /opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_BY

Rscript /opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/R/buildCRP.R in=/opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_BY/eqClass.txt out=/opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_BY/X_matrix.RData H=0.025

/opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/bin/TxIndexer -t /opt/working/projects/prj_054_1000seq_RNAseq/data/synII_new.transcripts.fa \
																			 -o /opt/working/projects/prj_054_1000seq_RNAseq/TxIndexer_idx_synII

Rscript /opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/R/genPolyesterSimulation.R /opt/working/projects/prj_054_1000seq_RNAseq/data/synII_new.transcripts.fa /opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_synII

/opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/bin/GenTC -i /opt/working/projects/prj_054_1000seq_RNAseq/TxIndexer_idx_synII \
																		 -l IU -1 /opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_synII/sample_01_1.fasta \
																		 -2 /opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_synII/sample_01_2.fasta \
																		 -o /opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_synII

Rscript /opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/R/buildCRP.R in=/opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_synII/eqClass.txt out=/opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_synII/X_matrix.RData H=0.025

/opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/bin/TxIndexer -t /opt/working/projects/prj_054_1000seq_RNAseq/data/BY4741.transcripts.fa \
																			 -o /opt/working/projects/prj_054_1000seq_RNAseq/TxIndexer_idx

Rscript /opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/R/genPolyesterSimulation.R /opt/working/projects/prj_054_1000seq_RNAseq/data/BY4741.transcripts.fa /opt/working/projects/prj_054_1000seq_RNAseq/design_matrix

/opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/bin/GenTC -i /opt/working/projects/prj_054_1000seq_RNAseq/TxIndexer_idx \
																		 -l IU -1 /opt/working/projects/prj_054_1000seq_RNAseq/design_matrix/sample_01_1.fasta \
																		 -2 /opt/working/projects/prj_054_1000seq_RNAseq/design_matrix/sample_01_2.fasta \
																		 -o /opt/working/projects/prj_054_1000seq_RNAseq/design_matrix

Rscript /opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/R/buildCRP.R in=/opt/working/projects/prj_054_1000seq_RNAseq/design_matrix/eqClass.txt out=/opt/working/projects/prj_054_1000seq_RNAseq/design_matrix/X_matrix.RData H=0.025


## creat eqc table
Rscript /opt/working/projects/prj_054_1000seq_RNAseq/eqc_code.R
bash /opt/working/projects/prj_054_1000seq_RNAseq/eqc.sh

## creat Ycount matrix
Rscript /opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/R/Create_count_matrix.R workdir=/opt/working/projects/prj_054_1000seq_RNAseq/eqc_BY design.matrix='/opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_BY/X_matrix.RData' core=8
Rscript /opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/R/Create_count_matrix.R workdir=/opt/working/projects/prj_054_1000seq_RNAseq/eqc_synII design.matrix='/opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_synII/X_matrix.RData' core=8
Rscript /opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/R/Create_count_matrix.R workdir=/opt/working/projects/prj_054_1000seq_RNAseq/eqc design.matrix='/opt/working/projects/prj_054_1000seq_RNAseq/design_matrix/X_matrix.RData' core=40

## AEM update Xmatrix
Rscript /opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/R/AEM_update_X_beta.R workdir=/opt/working/projects/prj_054_1000seq_RNAseq/eqc_BY core=8 design.matrix='/opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_BY/X_matrix.RData' isoform.out=XAEM_isoform_expression.RData paralog.out=XAEM_paralog_expression.RData merge.paralogs=FALSE isoform.method=average remove.ycount=TRUE
Rscript /opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/R/AEM_update_X_beta.R workdir=/opt/working/projects/prj_054_1000seq_RNAseq/eqc_synII core=8 design.matrix='/opt/working/projects/prj_054_1000seq_RNAseq/design_matrix_synII/X_matrix.RData' isoform.out=XAEM_isoform_expression.RData paralog.out=XAEM_paralog_expression.RData merge.paralogs=FALSE isoform.method=average remove.ycount=TRUE
Rscript /opt/working/projects/prj_054_1000seq_RNAseq/XAEM-binary-0.1.1/R/AEM_update_X_beta.R workdir=/opt/working/projects/prj_054_1000seq_RNAseq/eqc core=40 design.matrix='/opt/working/projects/prj_054_1000seq_RNAseq/design_matrix/X_matrix.RData' isoform.out=XAEM_isoform_expression.RData paralog.out=XAEM_paralog_expression.RData merge.paralogs=FALSE isoform.method=average remove.ycount=TRUE
