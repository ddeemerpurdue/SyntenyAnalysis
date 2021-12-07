'''
'''
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

myfile = 'rawdata/L-tmp-5000'


def main(file):
    uid = []
    starts = []
    ends = []
    distances = []
    lengths = []
    directions = []
    contigs = []
    count = 0
    prev_contig = ''
    prev_stop = ''
    less_zero = 0
    more_zero = 0
    zero = 0
    with open(file) as f:
        line = f.readline()
        line = f.readline()
        prev_contig = ''
        length = 0
        while line:
            gci, contig, start, stop, direction = line.split('\t')[0:5]
            start, stop = int(start), int(stop)
            if contig == prev_contig:
                contigs.append(contig)
                count += 1
                starts.append(start)
                ends.append(stop)
                uid.append(count)
                distance = start - prev_stop
                length = stop - start
                distances.append(distance)
                lengths.append(length)
                directions.append(direction)
                if distance < 0:
                    less_zero += 1
                elif distance == 0:
                    zero += 1
                elif distance > 0:
                    more_zero += 1
            else:
                pass
            prev_contig = contig
            prev_stop = stop
            line = f.readline()
    data = {'UID': uid, 'Start': starts, 'End': ends,
            'Distances': distances, 'Lengths': lengths,
            'Directions': directions, 'Contig': contigs}
    df = pd.DataFrame(data)
    print(
        f"Less than zero: {less_zero}\nZero: {zero}\nMore than zero: {more_zero}")
    return df


df = main(myfile)
print(df.head(5))

sns.scatterplot(x='Start', y='Distances', hue='Lengths', data=df)
# plt.show()
# sns.scatterplot(x='Lengths', y='Distances', hue='Directions', data=df,
#                 alpha=0.2)
# plt.show()

df_1500 = df[df['Lengths'] <= 1500]
# print(df_1500)
sns.scatterplot(x='Start', y='Distances', hue='Lengths', data=df_1500)
# plt.show()
df_large = df[df['Lengths'] > 1500]
# print(df_large)
sns.scatterplot(x='Start', y='Distances', hue='Lengths', data=df_large)
# plt.show()
print(len(df['Contig'].unique()))
graph = sns.FacetGrid(df, col='Contig', sharex=False)
graph.set(xlim=(0, 50000))
graph.map(plt.scatter, 'Start', 'Distances')
plt.show()
