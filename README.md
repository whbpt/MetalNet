# MetalNet

# pre install

pip install keras numpy graphviz requests
## more 
### for Ubuntu
sudo apt-get install graphviz
### for Centos
sudo yum install graphviz
### for Opensuse
sudo zypper install graphviz
# running script
#############
submit your target sequence in this webserver
http://gremlin.bakerlab.org/submit.php
after finish this calculation, please use the ID number as the script input;
for example, gremlin.bakerlab.org/submit.php&id=1535631343, the ID number is 1535631343
run the command
python predict.py 1535631343
