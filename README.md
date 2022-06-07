[![Netlify Status](https://api.netlify.com/api/v1/badges/2559146a-2438-4f28-a548-83786ff85ef2/deploy-status)](https://app.netlify.com/sites/web-lab-notes-hakyimlab/deploys)

## Lab notes in Rmarkdown

For this setup to work, you need to have blogdown and hugo installed on your computer. 
- `install.packages('blogdown') `
- `blogdown::install_hugo()`
- Open RStudio by double clicking `lab-notes.Rproj`
- Start a new analysis by adding a `New Post` from the addin option at the top of the `source` panel

<img src=https://github.com/hakyimlab/web-haky-personal-notes/blob/master/static/new-post-addin.png width="400x"> 

- Choose the title, author, tags
- In format select `R Markdown (.Rmd)` option
- git add, commit and push

Netlify is hosting the content  [here](https://lab-notes.hakyimlab.org)
