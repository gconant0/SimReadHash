#include "random.h"


Random_Gen::Random_Gen()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    
    seedin.open("seed.txt");
    if (seedin.fail()) {
        time (&rawtime);
        timeinfo = localtime(&rawtime);
        
        strftime(buffer,sizeof(buffer),"%A, %B %d, %Y %H:%M:%S",timeinfo);
        
        std::cout<<"No seed file found: seeding from time phrase "<<buffer<<std::endl;
        
        phrtsd(buffer, &seed1, &seed2);
    }
    else {
        seedin>>seed1>>seed2;
        seedin.close();
    }
    setall(seed1, seed2);
    
}

Random_Gen::Random_Gen(std::string seed_string)
{
    int i, stop, stop_seed;
    time_t rawtime;
    struct tm * timeinfo;
    char tbuffer[100], buffer[300];
    
    
    time (&rawtime);
    timeinfo = localtime(&rawtime);
        
    strftime(tbuffer,sizeof(tbuffer),"%A, %B %d, %Y %H:%M:%S",timeinfo);
    
    if (sizeof(buffer)<100) stop=sizeof(buffer);
    else stop=100;
    
    if (seed_string.length() >100) stop_seed=100;
    else stop_seed=seed_string.length();
    
    for(i=0; i<stop_seed; i++)
        buffer[i]=seed_string[i];
    
    for(i=0; i<stop; i++)
        buffer[i+stop_seed]=tbuffer[i];
    
    buffer[stop]='\0';
    
        
    phrtsd(buffer, &seed1, &seed2);
    
    setall(seed1, seed2);
    
}

Random_Gen::~Random_Gen()
{
    seedout.open("./seed.txt");
    getsd(&seed1, &seed2);
    seedout<<seed1<<"\t"<<seed2<<std::endl;
    seedout.close();
}
