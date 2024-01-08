import axios from 'axios';
import authHeader from './auth-header';

const API_URL = 'http://localhost:8080/api/v1/service/';

class UserService {
  getPublicContent() {
    const user = JSON.parse(localStorage.getItem('user'));
    return user ? axios.get(API_URL + 'user', { headers: authHeader(), params: {username:user.username} }) : {};
  }
  postCreateContainer() {
    const user = JSON.parse(localStorage.getItem('user'));
    return user ? axios.post(API_URL + 'create', {username:user.username}, { headers: authHeader()} ) : {};
  }
}

export default new UserService();
