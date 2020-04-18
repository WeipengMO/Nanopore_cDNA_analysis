import click


@click.command()
@click.option('--positive', required=True)
@click.option('--negative', required=True)
@click.option('--output', required=True)
def main(positive, negative, output):
    with open(positive, 'r') as pos, open(negative, 'r') as neg, open(output, 'w') as o:
        for posline, negline in zip(pos, neg):
            posparts=posline.rstrip().split("\t")
            negparts=posline.rstrip().split("\t")
            total=int(posparts[2])+int(negparts[2])
            o.write(f'{posparts[0]}\t{posparts[1]}\t{total}\n')


if __name__ == '__main__':
    main()