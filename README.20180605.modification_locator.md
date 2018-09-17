#	Modification Locator Processing

##	20180605

From last ubuntu image, create virtual machine.
*	name : ubuntu
* username : jake
* ssh public key : cat ~/.ssh/id_rsa.pub
* resource group : ubuntu


Using existing 20180604 with corrected evidence.txt



IP Address 40.121.39.184





```BASH
ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@40.121.39.184

sudo apt update
sudo apt full-upgrade
sudo apt autoremove
sudo reboot

ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@40.121.39.184

cd ~/syryu
git pull
make install

cd ~/modification_locator

modification_locator.rb --amino_acids STY --evidence evidence.txt --protein uniprot-organism+homo+sapiens.fasta > 20180605.modification_locator.txt &
```






###	UPLOAD TO AZURE STORAGE

Browse to portal.azure.com, Storage Account -> ryulab -> Access Keys to find a key.

Cleanup and upload data to Azure Storage and prep to save VM image ...

Remotely ...

```BASH
cd ~/modification_locator/
mkdir 20180605
mv MatchedModification*.txt ModificationNormalizerPeptides.txt ProteinModification.txt *.modification_locator.txt 20180605/
azcopy --verbose --source ~/modification_locator/20180605/ --destination https://ryulab.file.core.windows.net/ryulab/Modification%20Locator/20180605 --recursive --dest-key $( cat ~/dest-key )



sudo waagent -deprovision
```

Using the web portal GUI, save the image


