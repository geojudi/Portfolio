# Judith A. Sclafani - Portfolio

The follow collection contains samples of R scripts written for scientific projects. Examples have been selected to showcase my statistical and data processing capabilities.

## Table of Contents

## Building taxonomic abundance matrices
### **SKILLS:** dataframe and data file manipulation, data culling for quality and statistical robustness, function writing
This script is a part of a larger project to evaluate whether changes in abundance distributions through geologic time can be explained by neutral ecological theory. To do this, I wrote a program to transform a download from the Paleobiology Database into the format required by a previously developed maximum likehood estimation program written in PARI/GP. 

I needed to automate the analysis of thousands of fossil collections throughout geologic time. To do this, I developed a workflow, executed in bash, to cull the database download, split it by collection into multiple abundance matrices, save those matrices as .csv files, combine the abundance data with the additional information needed to run the PARI/GP program, excute the PARI/GP code, and combine the results into a new dataframe for additional analysis in R. 

Included here is the R script to perform the first part of this workflow, from culling the database download to writing the abundance matrices [PBDBtoAbundMatrices.r](/PBDBtoAbundMatrices.r)


## Ecological similiarity 


## Welcome to GitHub Pages

You can use the [editor on GitHub](https://github.com/geojudi/geojudi.github.io/edit/master/README.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/geojudi/geojudi.github.io/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
