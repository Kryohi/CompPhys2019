


int * currentTime()
{
    time_t now;
    struct tm *now_tm;
    static int currenttime[2];

    now = time(NULL);
    now_tm = localtime(&now);
    currenttime[0] = now_tm->tm_hour;
    currenttime[1] = now_tm->tm_min;
    
    return currenttime;
}


void make_directory(const char* name) 
{
    struct stat st = {0};
    
    #ifdef __linux__
       if (stat(name, &st) == -1) { mkdir(name, 0777); }
    #else
       _mkdir(name);
    #endif
}


void print_path()
{
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        printf("Current working dir: %s\n", cwd);
    } else {
        perror("getcwd() error");
    }
}
     
