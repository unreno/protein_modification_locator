#	Modification Locator Processing

##	20180601

From last ubuntu image, create virtual machine.
*	name : ubuntu
* username : jake
* ssh public key : cat ~/.ssh/id_rsa.pub
* resource group : ubuntu

Standard F16s_v2 (16 vcpus, 32 GB memory)

IP Address 23.101.142.194



```BASH
ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@23.101.142.194

sudo apt update
sudo apt full-upgrade
sudo apt autoremove
sudo reboot

ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@23.101.142.194

cd ~/syryu
git pull
make install

mkdir ~/modification_locator
cd ~/modification_locator

azcopy --verbose --source https://ryulab.file.core.windows.net/ryulab/Modification%20Locator/evidence.txt --destination evidence.txt --source-key $( cat ~/dest-key )
```


Locally, upload basic infrastructure ...

```BASH
rsync --archive --verbose --compress --rsh "ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress ~/github/unreno/syryu/modification_locator/PeptideMatch* jake@23.101.142.194:modification_locator/ 

rsync --archive --verbose --compress --rsh "ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress ~/github/unreno/syryu/modification_locator/uniprot-organism+homo+sapiens* jake@23.101.142.194:modification_locator/ 
```


```BASH
modification_locator.rb --amino_acids STY --evidence evidence.txt --protein uniprot-organism+homo+sapiens.fasta > 20180601.modification_locator.txt &
```


Wow is this slow. Pegs 1 CPU at 100% and using about 6GB memory. Just seems to sit there.

Quit. Need to figure out how to parallelize




###	UPLOAD TO AZURE STORAGE

Browse to portal.azure.com, Storage Account -> ryulab -> Access Keys to find a key.

Cleanup and upload data to Azure Storage and prep to save VM image ...

Remotely ...

```BASH
cd ~/modification_locator/
mkdir 20180601
mv MatchedModification*.txt ModificationNormalizerPeptides.txt ProteinModification.txt *.modification_locator.txt 20180601/
azcopy --verbose --source ~/modification_locator/20180601/ --destination https://ryulab.file.core.windows.net/ryulab/Modification%20Locator/20180601 --recursive --dest-key $( cat ~/dest-key )



sudo waagent -deprovision
```

Using the web portal GUI, save the image


