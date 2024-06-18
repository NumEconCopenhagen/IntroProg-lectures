from github import Github # pip install PyGithub


def load_repos(print_groups=False,check_year=False):
    # load personal access token
    with open(r"C:\Users\hms467\OneDrive - University of Copenhagen\Documents\Arbejde\Undervisning\IPNA_2024\Github\token.txt", mode = "r") as file:
        token = file.read()

    year = "2024"    
    class_name = "projects-" + year


    # a. access github through access token.
    gh = Github(token)
    org = gh.get_organization('NumEconCopenhagen')
    all_repos = org.get_repos()

    # b. locate all repositories in current class
    disregard_repos = ['lucas-og-anna']
    add_repos = ['Projects_2024_Emma-Anna-Oscar']

    current_class = [repo.name for repo in all_repos if ( (class_name in repo.name) & (repo.name not in disregard_repos) ) or (repo.name in add_repos)]

    if print_groups:
        # see this years' repos
        for r in current_class:
            print(r.removeprefix(class_name+"-"))


    if check_year:
        print(f'\nFollowing repos are not in the class but includes the year {year}:')

        for repo in all_repos:
            if (year in repo.name) & (repo.name not in current_class):
                print(repo.name)


    return class_name, all_repos, current_class, disregard_repos