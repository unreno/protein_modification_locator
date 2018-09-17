#	Modification Locator Processing

##	20180604

From last ubuntu image, create virtual machine.
*	name : ubuntu
* username : jake
* ssh public key : cat ~/.ssh/id_rsa.pub
* resource group : ubuntu





Using parallel ruby gem so gonna launch a VM with a bunch of cpus

Standard D64s v3 (64 vcpus, 256 GB memory)

(Could use E series for 432 GB memory but lets see how this goes)


Failed creation! Never had that happen on Amazon!

"message": "Allocation failed. We do not have sufficient capacity for the requested VM size in this region. Read more about improving likelihood of allocation success at http://aka.ms/allocation-guidance"

Trying Standard E64s v3 (64 vcpus, 432 GB memory)

Failed!

Trying Standard E32s v3 (32 vcpus, 256 GB memory)

Failed!

Trying Standard E64s v3 (64 vcpus, 432 GB memory) again

I deleted some empty resource groups just in case that was part of the problem?
Worked this time???



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

modification_locator.rb --amino_acids STY --evidence evidence.txt --protein uniprot-organism+homo+sapiens.fasta > 20180604.modification_locator.txt &
```






###	UPLOAD TO AZURE STORAGE

Browse to portal.azure.com, Storage Account -> ryulab -> Access Keys to find a key.

Cleanup and upload data to Azure Storage and prep to save VM image ...

Remotely ...

```BASH
cd ~/modification_locator/
mkdir 20180604
mv MatchedModification*.txt ModificationNormalizerPeptides.txt ProteinModification.txt *.modification_locator.txt 20180604/
azcopy --verbose --source ~/modification_locator/20180604/ --destination https://ryulab.file.core.windows.net/ryulab/Modification%20Locator/20180604 --recursive --dest-key $( cat ~/dest-key )



sudo waagent -deprovision
```

Using the web portal GUI, save the image


