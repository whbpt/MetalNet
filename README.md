# MetalNet 

MetalNet is designed for detecting metal binding sites with coevolution data. And metal binding sites and metal type will be offered if predicted.

# installation
```bash
pip install keras numpy graphviz requests networkx
```
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
1. submit your target sequence in this webserver
**http://gremlin.bakerlab.org/submit.php**
------------------
2. after finish this calculation, you will get a link such as, 
**http://gremlin.bakerlab.org/submit.php&id=1535637622**
------------------
3. remember to click the button,"**Import to GREMLIN BETA**".The ID number is 1535631343 (for example), and please use the ID number as the script input;
------------------ 
- **4. run the python script**
```bash
python predict.py 1535637622
```
