# MetalNet

# installation

pip install keras numpy graphviz requests
## one more step for graphviz 
- **for Ubuntu**
```bash
sudo apt-get install graphviz
```
- **for CentOS**
```bash
sudo yum install graphviz
```
- **for OpenSuse**
```bash
sudo zypper install graphviz
```
# running script
submit your target sequence in this webserver\n
http://gremlin.bakerlab.org/submit.php \n
after finish this calculation, please use the ID number as the script input;\n
for example, gremlin.bakerlab.org/submit.php&id=1535631343, the ID number is 1535631343\n
run the command \n
```bash
python predict.py 1535631343
```
