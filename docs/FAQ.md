# FAQ

### How can I ensure that my files have not been corrupted while being transferred ?

You can compute a hash value for each file and send that information along with your files
so people at the other end can confirm the data integrity.

In a Unix/Linux enviroment you can use the md5sum command. Here is an example:

```sh
$ md5sum *.gz > checksums.txt
$ cat checksums.txt
baa933d48c7dadb7a2ee46a67621e0f2  18277.vcf.anno.gz
d6bbba1fd1895979669abc56d2239d46  19466.vcf.anno.gz
84a8a85abe4f85e335143a3706d8a071  23138.vcf.anno.gz
06e9e55d0c41a3d80b103f7069ace4ba  27347.vcf.anno.gz
90c780176da6bb93686723b43300dc14  30111.vcf.anno.gz
c70aed11487f92ad795bbac2bf79747a  30114.vcf.anno.gz
```

You would transfer the checksums.txt file along with the rest of the files. Then, at the other
end, the collaborator would compute the hash values for the files and compare them with the
ones you sent:


```sh
$ md5sum -c checksums.txt
18277.vcf.anno.gz: OK
19466.vcf.anno.gz: OK
23138.vcf.anno.gz: OK
27347.vcf.anno.gz: OK
30111.vcf.anno.gz: OK
30114.vcf.anno.gz: OK
```
