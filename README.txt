docker build -t quay.io/benbfly/epiread_tools .
docker login https://quay.io/
docker push quay.io/benbfly/epiread_tools:latest

docker run -it -v '/mnt/bentest:/bentest' -v '/mnt/input:/input' -v '/home/azureuser/projects/nf-core:/nf-core' quay.io/benbfly/epiread_tools:latest /bin/bash

# Tools for modbam_parser
~/projects/nf-core/htslib/test/base_mods$ ../test_mod /mnt/bentest/testdata/nanoseq/modbam2bed_berman/test-fw-rev.sam

# Mm/Ml tags are stranded, so if the read is reverse they start at the end of the sequence relative to the BAM
python3 ./epiread_tools/epireadToBedgraph.py -j /nf-core/epiread-tools/tests/modbam_test.json



# Deconvolution
cd ../deconvolution_models

# BAM
# NOTE: IF YOU RUN /opt/venv/bin/deconvolution, it will still get the /opt/venv libraries
/nf-core/epiread-tools/deconvolution_models/deconvolution --model celfie-ish -a /input/Data/Methylation/Unterman_2023/hg38/beta_atlas.hg38.tsv \
-i /input/Data/Methylation/Unterman_2023/hg38/U250.hg38.tsv -b \
-m /nf-core2/testresults_good/mergedModBamWithDupsMarked/INV_0002_fake.bam --epiformat modbam \
--num_iterations 1 --stop_criterion 0.0000001 \
--random_restarts 1 --outfile ~/deconvolution_demo_output.txt \
--cpg_coordinates /input/Data/Genomes/hg38_noalt/hg38.analysisSet.CGmotifs.modkit.merged.1based.bed.gz \

# DEMO
# NOTE: IF YOU RUN /opt/venv/bin/deconvolution, it will still get the /opt/venv libraries
/nf-core/epiread-tools/deconvolution_models/deconvolution --model celfie-ish -a demo/beta_atlas.txt \
-i demo/U250.tsv -b \
-m demo/mixture.epiread.gz --epiformat old_epiread_A \
--num_iterations 1 --stop_criterion 0.0000001 \
--random_restarts 1 --outfile ~/deconvolution_demo_output.txt \
--cpg_coordinates demo/sample_cpg_file.bed.gz \
#-b \
#--cpg_coordinates demo/hg19.CpG.bed.sorted.gz \

# Working deconvolution
# NOTE: IF YOU RUN /opt/venv/bin/deconvolution, it will still get the /opt/venv libraries
/nf-core/epiread-tools/deconvolution_models/deconvolution \
--model celfie-ish \
-a /input/Data/Methylation/Unterman_2023/hg38/beta_atlas.hg38.tsv \
-b -i /input/Data/Methylation/Unterman_2023/hg38/U250.hg38.tsv \
-m /bentest/modbams-Katsman-2022/ISPRO.bc11.meg242.remora1.edgefix.mod_mappings.sorted.hg38.bam \
 --epiformat modbam --num_iterations 1000 --stop_criterion 0.0000001 --random_restarts 1 \
 --outfile ~/deconvolution_demo_4output.txt \
 --cpg_coordinates /input/Data/Genomes/hg38_noalt/hg38.analysisSet.CGmotifs.modkit.merged.0based.bed.gz

 # UXM demo
 # NOTE: IF YOU RUN /opt/venv/bin/deconvolution, it will still get the /opt/venv libraries
/nf-core/epiread-tools/deconvolution_models/deconvolution --model uxm -i demo/U250.tsv -b -m demo/mixture.epiread.gz --epiformat old_epiread_A \
 --random_restarts 1 --outfile ~/deconvolution_demo_output_uxm.txt --cpg_coordinates demo/sample_cpg_file.bed.gz \
 --min_length 2 --u_threshold 0.25 --percent_u demo/U_atlas.txt
 
# Epistate demo (randome restarts, num_iterations, stop_criterion reqiured) NOTE: --num_iterations 10000 --stop_criterion 0.0000001 takes quite a long time
# NOTE: IF YOU RUN /opt/venv/bin/deconvolution, it will still get the /opt/venv libraries
/nf-core/epiread-tools/deconvolution_models/deconvolution --model epistate -i demo/U250.tsv -b -m demo/mixture.epiread.gz --epiformat old_epiread_A \
--outfile ~/deconvolution_demo_output_epistate.txt --cpg_coordinates demo/sample_cpg_file.bed.gz \
--lambdas demo/lambdas.txt --thetas demo/thetas.txt --random_restarts 1 --num_iterations 100 --stop_criterion 0.0001

# UXM nanopore (NOT CURRENTLY WORKING)
deconvolution --model uxm -b -i /input/Data/Methylation/Unterman_2023/hg38/U250.hg38.tsv \
-m /bentest/modbams-Katsman-2022/ISPRO.bc11.meg242.remora1.edgefix.mod_mappings.sorted.hg38.bam \
 --epiformat modbam --random_restarts 1  --outfile ~/deconvolution_modbam_uxm_restart1_minlen25.txt \
 --cpg_coordinates /input/Data/Genomes/hg38_noalt/hg38.analysisSet.CGmotifs.modkit.merged.0based.bed.gz \
 --min_length 2 --u_threshold 0.25 --percent_u /input/Data/Methylation/Unterman_2023/hg38/U_atlas.hg38.tsv

# Epistate nanopore (NOT CURRENTLY WORKING MAYBE HG38 problem)
# (not sure iterations. --num_iterations 10000 --stop_criterion 0.0000001 takes quite a long time)
# NOTE: IF YOU RUN /opt/venv/bin/deconvolution, it will still get the /opt/venv libraries
/nf-core/epiread-tools/deconvolution_models/deconvolution --model epistate -b -i /input/Data/Methylation/Unterman_2023/hg38/U250.hg38.tsv \
 -m /bentest/aug4_HL_stratified_merge_012/mergedModBamRemovedSupsAndDups/ONT-R005-B8-barcode09.nosupsnodups.bam \
 --epiformat modbam --random_restarts 1  --outfile ~/deconvolution_modbam_epistate_restart1_minlen25.txt \
 --cpg_coordinates /input/Data/Genomes/hg38_noalt/hg38.analysisSet.CGmotifs.modkit.merged.0based.bed.gz \
 --lambdas /input/Data/Methylation/Unterman_2023/hg38/lambdas.U250.hg38.tsv \
 --thetas /input/Data/Methylation/Unterman_2023/hg38/thetas.U250.hg38.tsv \
 --random_restarts 1 --num_iterations 100 --stop_criterion 0.0001


# Is there a problem with this bam?  Something incompatible in the old bam format?
# -m /bentest/modbams-Katsman-2022/ISPRO.bc11.meg242.remora1.edgefix.mod_mappings.sorted.hg38.bam \
