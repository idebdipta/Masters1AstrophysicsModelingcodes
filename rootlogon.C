{
  gSystem->Load("libTree.so");
  gSystem->Load("libHist.so");
  gSystem->Load("libGpad.so");
  gSystem->Load("/usr/lib/x86_64-linux-gnu/libgsl.so");
  gSystem->SetMakeSharedLib("cd $BuildDir ; g++ -c $Opt -pipe -Wall -W -Woverloaded-virtual -fPIC -O3 -g -Iinclude -pthread $IncludePath $SourceFiles ;  g++ $ObjectFiles -shared -Wl,-soname,$LibName.so -O  -g -o $SharedLib");
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
}



