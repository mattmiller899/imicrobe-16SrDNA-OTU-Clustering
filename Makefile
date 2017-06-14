APP = imicrobe-16SrDNA-OTU-Clustering
VERSION = 0.0.1
EMAIL = jklynch@email.arizona.edu

clean:
	find . \( -name \*.conf -o -name \*.out -o -name \*.log -o -name \*.param -o -name launcher_jobfile_\* \) -exec rm {} \;

files-delete:
	files-delete $(CYVERSEUSERNAME)/applications/$(APP)-$(VERSION)

files-upload:
	files-upload -F stampede/ $(CYVERSEUSERNAME)/applications/$(APP)-$(VERSION)

apps-addupdate:
	apps-addupdate -F stampede/app.json

deploy-app: clean files-delete files-upload apps-addupdate

test:
	sbatch test.sh

jobs-submit:
	jobs-submit -F stampede/job.json

container:
	rm -f singularity/$(APP).img
	sudo singularity create --size 2048 singularity/$(APP).img
	sudo singularity bootstrap singularity/$(APP).img singularity/$(APP).def

iput-container:
	iput -K singularity/$(APP).img

iget-container:
	iget -K $(APP).img
	mv $(APP).img stampede/

