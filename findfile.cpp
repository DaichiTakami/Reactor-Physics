#include<iostream>
#include<filesystem>
namespace fs = std::filesystem;

int main(int argc, char** argv){

	if (argc == 1)return -1;
	fs::path path1(argv[1]);

	if(fs::is_directory(path1)){
		std::cout<<"Directory:"<<fs::absolute(path1).string()<<std::endl;
		std::cout<<"Files:"<<std::endl;
		auto dir_it = fs::directory_iterator(path1);
		for(auto &p : dir_it){
			std::cout<<p.path().string()<<std::endl;
		}
	}

	auto path2 = path1/"new_dir";
	fs::create_directory(path2);
	
}