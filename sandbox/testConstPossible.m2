


createFoo = ()->
(
  foo := new MutableHashTable;
  fkt := ()->null;
  foo.bar = ()->     (    return fkt();  );
  foo.setBar = (m)-> (    fkt = copy m;  );
  return new HashTable from foo;
)


----------------------------
foo = createFoo();
trick = ()->5;
bar = ()->trick();
foo.setBar(bar)
foo.bar()
trick = ()->0;
foo.bar()  --tricked!


