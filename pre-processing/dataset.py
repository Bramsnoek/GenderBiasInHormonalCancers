import requests

temp = []
type_ages = {"Cervix": [18, 30], "Ovary": [18, 40], "Pancreas": [18, 50], "Prostate": [40, 50], "Testis": [18, 25],
             "Thyroid": [18, 30]
             }


def get_cases(location):
    """
    Deze functie haalt samples en bijbehorende metadata op.
    :param location: Kankerssoort
    :return: Lijst met samples
    """
    print("Starting", location)
    fields = [
        "cases.demographic.race",
        "file_name",
        "cases.samples.sample_type",
        "cases.disease_type",
        "cases.project.primary_site",
        "cases.demographic.gender",
        "cases.files.date_release",
        "cases.diagnoses.age_at_diagnosis"

    ]
    fields = ",".join(fields)
    files_endpt = "https://api.gdc.cancer.gov/files"

    # This set of filters is nested under an 'and' operator.
    # Leeftijd, uploaddatum, ethniticity
    filters = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.primary_site",
                    "value": [f"{location}"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.experimental_strategy",
                    "value": ["RNA-Seq"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.type",
                    "value": ["gene_expression"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_format",
                    "value": ["TXT"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "cases.samples.sample_type",
                    "value": ["Primary Tumor"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "cases.demographic.gender",
                    "value": ["Male", "Female"]
                }
            }
        ]
    }

    # A POST is used, so the filter parameters can be passed directly as a Dict object.
    params = {
        "filters": filters,
        "fields": fields,
        "format": "TSV",
        "size": "200000"
    }

    response = requests.post(files_endpt, headers={"Content-Type": "application/json"}, json=params)
    cases = []
    data = response.content.decode()
    male = 0
    female = 0
    try:
        for case in data.split("\n"):
            split_data = case.split("\t")
            if split_data[6].__contains__("htseq.counts"):
                try:
                    split_data[2] = round(int(split_data[2]) / 365, 1)
                    split_data[1] = split_data[1].replace(' ', '_')
                    if type_ages.get(location)[0] < int(split_data[2]) <= type_ages.get(location)[1] and split_data[1] \
                            == 'white':
                        cases.append(split_data)
                        if split_data[0] == 'male':
                            male += 1
                        else:
                            female += 1
                except ValueError:
                    pass
    except IndexError:
        pass
    print(
        f"{location} with sample aged {type_ages.get(location)[0]}-{type_ages.get(location)[1]}: \nMale:\t{male}\nFemale:\t{female}\n")
    return cases


def write_tsv(case):
    """
    Deze functie schrijft metadata weg in een TSV format.
    :param case: Sample met metadata
    """
    with open("metadata.tsv", "a") as tsv_file:
        file_case = case[6].strip(".gz")
        age_cat = f"{type_ages.get(case[4])[0]}-{type_ages.get(case[4])[1]}"
        tsv_file.write(f"{file_case}\t{case[4]}\t{case[0]}\t{case[2]}\t{age_cat}\n")


def write_data(case_by_id, loc):
    """
    Verkrijgt de HTSeq counts met sample ID's
    :param case_by_id: Een verkregen lijst van de get_cases functie. Hierin staat sampledata
    :param loc: Kankersoort
    """
    for fetch_data in case_by_id:
        data_endpt = "https://api.gdc.cancer.gov/data/{}".format(fetch_data[7].rstrip())
        response = requests.get(data_endpt, headers={"Content-Type": "application/json"})

        try:
            # Bepaal leeftijd, niet mijn finest work of code :')
            age = 'Not reported'
            # Bestanden worden geschreven naar folder: dataset/kankersoort/geslacht/leeftijdscategorie/ethniciteit/
            with open(f"dataset/{loc}/{fetch_data[6].rstrip()}",
                      "wb") as output_file:
                print(
                    f"Downloading case: dataset/{loc}/{fetch_data[6].rstrip()}"
                )
                output_file.write(response.content)
                write_tsv(fetch_data)
        except Exception as e:
            print(e)
            print(f"{fetch_data[6]} could not be downloaded, skipping...")


def main():
    with open("metadata.tsv", "w") as file:
        file.write("Case_id\tPrimary_site\tGender\tAge\tAge category\n")
        file.close()
    for loc in type_ages.keys():
        data = get_cases(loc)
        write_data(data, loc)


main()
