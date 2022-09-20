
# activate conda env
source activate yevo_pipeline_env

# run pipeline
snakemake --snakefile workflow/Snakefile.py --configfile config/config.yml --cores 4 --use-conda







