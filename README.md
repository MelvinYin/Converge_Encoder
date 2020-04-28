# converge_encoder

Supports `converge` by parsing and converting `proteome.fasta`, 
`seed_seqs.fasta`, and `BLOSUM62` into int representations. 

Output as binary archives (`proteome_binary`, `seed_seq_binary`, 
`blosum_binary`) in cereal format. 

C++17, cereal v1.2.2, cmake. 

Requirements: <br>
brew install build-essential cmake openmpi

Docker link: <br>
todo

To run: <br>
Local
* cd to project directory.
* `cmake ./`
* `make`
* Move `proteome.fasta`, `seed_seqs.fasta` into `input`.
* `./converge_encoder`
* Take `proteome_binary`, `seed_seq_binary`, `blosum_binary` from `output`. <br>

Docker
* Move `proteome.fasta`, `seed_seqs.fasta` into `input`.
* docker run -i docker_image_name
* [local] docker cp ./converge_encoder 66de272fe9c9:/home
* [docker] cd /converge_encoder
* [docker] cmake ./
* [docker] make
* [docker] ./converge_encoder
* [local] docker cp 66de272fe9c9:/home/converge_encoder/output ./

To save:
* [local] docker commit 66de272fe9c9 docker_image_name

To run:
* [local] docker run -i docker_image_name

Expected runtime 2 seconds. 
