# AiiDA Vasp Plugin

## Install and usage:
### Install AiiDA
1. Download from [www.aiida.net/?page_id=264](www.aiida.net -> Download) (only EPFL version supports ssh)
2. install and setup -> [http://aiida-core.readthedocs.org/en/stable/](aiida's documentation)

### Install this plugin
From the aiida_vasp folder use:

```bash
$ python install.py <path to aiida distribution folder>
# or type
$ python install.py -h # for more info and install options
```

If you are going to try and use the current version for production, it's recommended to

```bash
$ python install.py --copy <aiida folder> vasp_XX # where XX is some ID
```

And use exclusively the vasp_XX plugins for production, so that changes in later versions 
don't make your work harder than necessary.

To test wether the installation was successful use

```bash
$ verdi calculation plugins 

# example output:

## Pass as a further parameter one (or more) plugin names
## to get more details on a given plugin.
* codtools.cifcellcontents
* codtools.cifcodcheck
* codtools.cifcodnumbers
* codtools.ciffilter
* codtools.cifsplitprimitive
* quantumespresso.cp
* quantumespresso.pw
* quantumespresso.pwimmigrant
* simpleplugins.templatereplacer
* vasp.asevasp
* vasp.vasp
```
You should see vasp.asevasp and vasp.vasp in the list

### Configure the code
See [http://aiida-core.readthedocs.org/en/stable/setup/computerandcodes.html#computer-setup-and-configuration](aiida docs)
on how to set up computers and codes. Note that you need at least one computer configured and a VASP executable on it
in order to use this plugin.

### Using this plugin
in aiida_vasp/tests there are some example python scripts which submit basic vasp runs as well as a basic bandstructure workflow. 
Take a look at the aiida docs in order to understand what exactly they do.
