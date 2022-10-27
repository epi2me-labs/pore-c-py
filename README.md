## Contributing

Prepare your environment

```shell
pipx install phatch       # to manage dependencies & environments
pipx install pre-commit   # clean code
pipx install commitizen   # clean commit messages
pipx install tox          # test against multiple versions of python
```

Initial setup

```shell
pre-commit install
pre-commit install --hook-type commit-msg --hook-type pre-push
hatch run pore-c2 --help        # create a virtualenv for testing
```

Developer Workflow

```shell

# 1. Create a branch
git checkout -b <branch-name>
# 2. Edit some files
# 3. Run some tests
hatch run test    # run pytest tests
# 4. Get ready to commit your work
git add <edited files>  # stage your files
git commit      # this will run the precommit hooks
# 5. Fix any things flagged by pre-commit that haven't been autofixed
# 6. Once things have passed use commitizen to create a nice commit messaged
cz commit
# 7. Squash and merge
git rebase -i main
git checkout main
git merge <branch-name>
```
