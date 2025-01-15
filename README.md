![Rib](https://github.com/user-attachments/assets/61a4005a-3d11-4fe3-9e28-b8efb75c747b)

## Authors
* **Milo Scola**(<miloscola@berkeley.edu>) - *University of California Berkeley*

* **Professor Christopher Anderson** - *University of California Berkeley*

## Function
This program selects the optimal ribosome binding site from a list of genetic parts for a gene to minimize hairpins, ensuring efficent translation. 
This program is designed to be compatable with Professor Christopher Anderson's SYNISTER database and can be found under composition/RBSchooser2.

## Usage
1) add a list of RBS parts in your lab to composition/data/rbs_options in TSV format: name, description, RBS, CDS
2) input the construct you are adding a RBS to to RBSchooser2 and run the program. The program will return the RBS which minimizes the number of hairpins. 


