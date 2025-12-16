
### Branch Transcriptomic (annotate_cells): api call(optional; api_hca_userinp.py) --> normalization(scanpy_norm.py) --> annotation(annotate_celltypes.py)

*This is (for the time being) seperated as a branch, because some processed matrix data from HCA contain annotation already (hopefully). In that case it'd be preferable to choose those*

`scanpy_norm.py` :-

The idea of this code is to prepare single-cell data for downstream analysis like the standard workflow. 

Normalizes each cell to a fixed count, then applies log1p transformation

THIS IS NEEDED FOR CELLTYPIST-ANNOTATION

`annotate_celltypes.py` :-

This is to have a column in the loom/h5ad file (AnnData format), storing the cell subtypes eg. -> 
![Description][/home/biodata/aman/scripts/Pasted image 20251202131332.png]

in a layer called 'cell_subtype' in `obs` layer (changeable)

![Description][/home/biodata/aman/scripts/Pasted image 20251202131416.png]

So that in the next steps, pooling as suggested by [*Gustafsson et al.*](https://doi.org/10.1073/pnas.2217868120) can be done according to the cell types

---

### Branch Transcriptomic (metabolic):  api call(optional; api_hca_userinp.py) --> normalization(norm_pooling.py) --> gene symbol conversion (genetoensembl.py) --> contextualization algorithm (gimme.py | tinit.py | fastcore.py)

`norm_pooling.py` :-

Ideas is to aggregate cells by type. Pools or sums raw counts across cells (either all cells or by cell type; psuedo-bulk), then converts to counts per million; *Gustafsson et al.* -

![Description][/home/biodata/aman/scripts/SCR-20251202-mpyd.png]

`genetoensembl.py` :-

This code will read a gtf file and maps the expression data to ensembles gene ids. The version number is stripped out. To be compatible with Human-GEM model

`gimme.py` :-

`tinit.py` :-

`fastcore.py` :-

`pymCADRE.py` :-

---