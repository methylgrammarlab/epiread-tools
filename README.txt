docker build -t epiread_tools .
docker run -it -v '/mnt/bentest:/bentest' -v '/mnt/input:/input' -v '/home/azureuser/projects/nf-core:/nf-core' epiread_tools