O = "./R/objects"

all: data peptides proteins

data: ./output/load.Rout ./output/clean.Rout
peptides: ./output/get_peptides2.Rout
proteins: ./output/get_peptides2.Rout ./output/peptides_correlations.Rout ./output/get_proteins.Rout
analysis1: ./output/get_peptides2.Rout ./output/peptides_correlations.Rout ./output/get_proteins.Rout ./output/analysis1.Rout

#loading data
./output/load.Rout: ./R/boot.R ./R/functions.R ./R/load.R ./data/*
	R CMD BATCH --no-save --no-restore ./R/load.R ./output/load.Rout

#cleaning data
./output/clean.Rout: ./R/clean.R ./R/objects/*_load_.RData
	R CMD BATCH --no-save --no-restore ./R/clean.R ./output/clean.Rout

#getting peptides filtered based on Qvalue from spectronaut also adjusting for batch effects

./output/get_peptides2.Rout: ./R/src/get_peptides2.R ./R/objects/peptides.data._clean_.RData ./R/objects/exp_metadata._clean_.RData
	R CMD BATCH --no-save --no-restore ./R/src/get_peptides2.R ./output/get_peptides2.Rout

#getting correlations among peptides within protein
./output/peptides_correlations.Rout: ./R/src/peptides_correlations.R ./R/objects/peptides.peak_sums.trimmed.RData ./R/objects/protein_annotations._load_.RData
	 R CMD BATCH --no-save --no-restore ./R/src/peptides_correlations.R ./output/peptides_correlations.Rout

#calculating protein fold-changes
./output/get_proteins.Rout: ./R/src/get_proteins.R ./R/objects/exp_metadata._clean_.RData ./R/objects/peptides.matrix.RData ./R/objects/peptides.matrix.combat.RData ./R/objects/peptides.cor.stats.top.RData ./R/objects/protein_annotations._load_.RData
	R CMD BATCH --no-save --no-restore ./R/src/get_proteins.R ./output/get_proteins.Rout

./R/objects/proteins.matrix.combat.quant.FC.RData: ./output/get_proteins.Rout

./R/objects/proteins.matrix.combat.quant.RData: ./output/get_proteins.Rout

./output/analysis1.Rout: ./R/src/analysis1.R ./R/objects/proteins.matrix.combat.quant.FC.RData ./R/objects/proteins.matrix.combat.quant.RData ./R/objects/exp_metadata._clean_.RData
	R CMD BATCH --no-save --no-restore ./R/src/analysis1.R ./output/analysis1.Rout


#./output/get_peptides.Rout: ./R/src/get_peptides.R ./R/objects/peptides.data.RData ./R/objects/experiment.map._load_.RData ./R/objects/dates_map._load_.RData \
#                            ./R/objects/sample.map._load_.RData ./R/objects/protein_annotations._load_.RData
#	R CMD BATCH --no-save --no-restore ./R/src/get_peptides.R ./output/get_peptides.Rout

#calculating correlations among peptides 
./output/peptides_correlations.Rout: ./R/src/peptides_correlations.R ./R/objects/peptides.peak_sums.trimmed.RData ./R/objects/protein_annotations._load_.RData
	R CMD BATCH --no-save --no-restore ./R/src/peptides_correlations.R ./output/peptides_correlations.Rout

#adjusting for batch effects on protein level
#./output/batch_effects.Rout: ./R/src/batch_effects.R ./R/objects/peptides.peak_sums.trimmed.RData ./R/objects/protein_annotations._load_.RData \
#                             ./R/objects/peptides.cor.stats.top.RData ./R/objects/sample_map._load_.RData
#	R CMD BATCH --no-save --no-restore ./R/src/batch_effects.R ./output/batch_effects.Rout



clean:
	rm ./output/*
