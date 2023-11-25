from selenium import webdriver
from selenium.webdriver.common.by import By
import time
import sys
import os


if __name__ == '__main__':

    # locus_name = 'D10S1248'
    locus_name = sys.argv[1]

    if os.path.exists(f'results/{locus_name}.genotype.txt'):
        print('Exist')
        sys.exit(0)


    # driver = webdriver.Chrome()
    driver = webdriver.Edge()
    try:
        driver.get(f'https://strbase.nist.gov/Human/VariantalleleTable/{locus_name}')
    except:
        print( f'{locus_name} can not reach' )

    max_attempts = 5
    attempts = 0

    while attempts < max_attempts:
        try:
            time.sleep(300)
            rows = driver.find_element(By.XPATH, '//*[@id="vaTable"]')
            break
        except:
            attempts += 1
            driver.get(f'https://strbase.nist.gov/Human/VariantalleleTable/{locus_name}')

    try:

        output = []
        for row in rows.text.split('\n'):
            if row[0].isdigit():
                col_lst = row.split(' ')
                output.append( locus_name + ',' + col_lst[0] + '\n')
        
        file = open(f'results/{locus_name}.genotype.txt', 'w')
        for o in output:
            file.write( o )
        file.close()

    except:
        print(f'{locus_name} page not open')

    driver.close()