# my-quarto-website

A friendly description of your project goes here.

## Usage

[Create a new repository](https://github.com/organizations/gandallab/repositories/new)
under the gandallab org and select this repo as a template. Be sure to enable
"Include all branches" (otherwise you will have to follow the instructions
[here](https://quarto.org/docs/publishing/github-pages.html#publish-command)
to create a gh-pages branch).

You can render individual notebooks by clicking the Render button in RStudio.

To render the whole website, switch to the terminal and run `quarto render`.

If your website renders successfully, you are ready to deploy to GitHub Pages.
Commit any changes to the main branch, then run `quarto publish gh-pages`.
This will update the gh-pages branch with the contents of _site.
