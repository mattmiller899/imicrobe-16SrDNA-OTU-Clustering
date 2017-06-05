APP = imicrobe-JUITS16S
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
	rm -f singularity/imicrobe-JUITS16S.img
	sudo singularity create --size 2048 singularity/imicrobe-JUITS16S.img
	sudo singularity bootstrap singularity/imicrobe-JUITS16S.img singularity/imicrobe-JUITS16S.def

iput-container:
	irm imicrobe-JUITS16S.img
	iput -K singularity/imicrobe-JUITS16S.img

iget-container:
	iget -K imicrobe-JUITS16S.img
	mv imicrobe-JUITS16S.img stampede/

