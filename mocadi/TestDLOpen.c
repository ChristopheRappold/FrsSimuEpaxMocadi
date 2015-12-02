#include <stdlib.h>
#include <dlfcn.h> 
#include <stdio.h>

int main()
{
  
  char plugin [] = {"/u/crappold/frs/mocadi/libFuncLoadTree.so"};
  char func [] = {"FuncLoadTree"};
  void *handle = dlopen(plugin, RTLD_LAZY);
  if(!handle)
    {
      printf("error: the external library %s was not found\n",plugin);
      printf("%s\n",dlerror());
      exit(1);
    }
  printf("any error? %s\n",dlerror());
  void* symbol = dlsym(handle,func);
  printf("symbol loading any error? %s\n",dlerror());

  return 0;
}
