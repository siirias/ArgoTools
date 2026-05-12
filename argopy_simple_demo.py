from argopy import DataFetcher as Argo
import matplotlib.pyplot as plt

poiju_wmo = '4903797'
profiili_n = 15
# Määritetään haettava WMO ja hakutapa (erddap on yleisin/vakio)
fetcher = Argo().float(poiju_wmo)

# Haetaan data xarray-muodossa
ds = fetcher.to_xarray()

# 2. Luodaan kuvaaja
plt.figure(figsize=(8, 6))

# Ryhmitellään data sykleittäin ja piirretään esim. 5 ensimmäistä
for cycle, data in ds.groupby('CYCLE_NUMBER'):
    if cycle > profiili_n: 
        break # Lopetetaan profiili_n kuvaajan jälkeen.
    plt.plot(data.TEMP, data.PRES, label=f'Sykli {cycle}')

# 3. Muotoilu
plt.gca().invert_yaxis() # Pinta ylös
plt.xlabel('Lämpötila (°C)')
plt.ylabel('Paine (dbar) / Syvyys')
plt.title(f'Argo {poiju_wmo}')
plt.legend()
plt.grid(True)

plt.show()


# 1. Luodaan kuvaaja
plt.figure(figsize=(12, 6))

# 2. Piirretään pisteet: x=aika, y=paine, väri (c)=lämpötila
# s=1 tekee pisteistä pieniä, jotta ne muodostavat jatkuvan pinnan
sc = plt.scatter(ds.TIME, ds.PRES, c=ds.TEMP, s=2, cmap='icefire')

# 3. Muotoilu
plt.gca().invert_yaxis()  # Pinta ylhäällä
plt.colorbar(sc, label='Lämpötila (°C)') # Lisätään väriselitys
plt.ylabel('Syvyys (dbar)')
plt.title(f'Argo {poiju_wmo}')

plt.show()