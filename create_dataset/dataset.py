import requests

temp = []


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
    try:
        for case in data.split("\n"):
            split_data = case.split("\t")
            if split_data[6].__contains__("htseq.counts"):
                try:
                    split_data[2] = round(int(split_data[2]) / 365, 1)
                    split_data[1] = split_data[1].replace(' ', '_')
                    cases.append(split_data)
                except ValueError:
                    pass
    except IndexError:
        pass
    return cases


def write_tsv(case):
    """
    Deze functie schrijft metadata weg in een TSV format.
    :param case: Sample met metadata
    """
    with open("metadata.tsv", "a") as tsv_file:
        tsv_file.write(f"{case[7].rstrip()}\t{case[5]}\t{case[4]}\t{case[0]}\t{case[1]}\t{case[2]}\n")


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
            if 0 <= fetch_data[2] <= 20:
                age = '0-20'
            elif 20 < fetch_data[2] <= 40:
                age = '21-40'
            elif 40 < fetch_data[2] <= 60:
                age = '41-60'
            elif 60 < fetch_data[2]:
                age = '61+'
            # Bestanden worden geschreven naar folder: dataset/kankersoort/geslacht/leeftijdscategorie/ethniciteit/
            with open(f"dataset/{loc}/{fetch_data[0]}/{age}/{fetch_data[1]}/{fetch_data[7].rstrip()}.gz",
                      "wb") as output_file:
                print(
                    f"Downloading case: dataset/{loc}/{fetch_data[0]}/{age}/{fetch_data[1]}/{fetch_data[7].rstrip()}.gz"
                )
                output_file.write(response.content)
                write_tsv(fetch_data)
        except Exception as e:
            print(e)
            print(f"{fetch_data[7]} could not be downloaded, skipping...")


def main():
    with open("metadata.tsv", "w") as file:
        file.write("Case_id\tSample_type\tPrimary_site\tGender\tAge\tEthnicity\n")
        file.close()
    locations = ["Pancreas", "Cervix", "Ovary", "Prostate", "Testis", "Thyroid"]
    for loc in locations:
        data = get_cases(loc)
        write_data(data, loc)


main()
