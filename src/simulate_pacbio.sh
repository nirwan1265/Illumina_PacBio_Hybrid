# Clone the git repository
# > git clone https://github.com/yukiteruono/pbsim3.git
# > ./configure
# > make
# > make install
# > pbsim -h

pbsim --strategy wgs --method errhmm --errhmm data/ERRHMM-RSII.model --depth 10 --length-min 10000 --length-max 12000 --genome simulated_genome/simulate_genome_500kb.fa

# strategy = wgs
# method = errhmm (error model)
# genome = fasta file
# depth = coverage (max 20.0)
# length-min = ( min 100)
# length-max = (max 1000000)
# --accuracy-min = ( min 0.75)

# For all files in a folder:
#!/bin/zsh

# set the input and output directory paths
input_dir="reference"
output_dir="."

# loop through each file in the input directory
for file in "${input_dir}"/*
do
  # get the filename without the path or extension
  filename=$(basename "${file}")
  filename="${filename%.*}"

  # run the command with the current input file and output filename
  pbsim --strategy wgs --method errhmm --errhmm data/ERRHMM-RSII.model --depth 10 --length-min 10000 --length-max 12000 --genome "${file}" --prefix "${output_dir}/simulated_${filename}_"
done
