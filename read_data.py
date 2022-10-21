import pandas as pd

# Datos considerados del 22-03-2020 al 07-09-2022 (900 dias en total)

df = pd.read_excel("dados_full_brasil.xlsx")

initial_date = 25 # 22-03-2020
total_days = 900

with open("output/data.txt","w") as f:
    f.write("%i\n" % total_days)
    for i in range(initial_date,initial_date + total_days):
        x = df["new_deaths_smoothed_per_million"][i]

        if pd.isna(x) == True:
            f.write("%f\n" % 0.0)
        else:
            f.write("%f\n" % x)
